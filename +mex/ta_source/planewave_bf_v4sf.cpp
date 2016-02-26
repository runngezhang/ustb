#include "planewave_bf.hpp"
#include <mex.h>

#define USE_SSE2
#include "sse_mathfun.h"

#ifdef _HAS_SSE4_
    #include <smmintrin.h>
#endif

#include <cmath>
#include <algorithm>

#define DEF_CONST(a,b) const v4sf a = _mm_set1_ps(b);
#define DEF_CONSTI(a,b) const v4si a = _mm_set1_epi32(b);

typedef union 
{
  __m128 m128f;
  __m128i m128i;  
  float f[4];
  int i[4];
} vec4;

typedef v4sf (*APOD_FUNCTION_V4SF)(v4sf);
typedef v4sf (*XX_FUNCTION)(v4sf,v4sf);

DEF_CONST(CF4_0,        0.0f) 
DEF_CONST(CF4_1,        1.0f) 
DEF_CONST(CF4_1000,        1000.0f) 
DEF_CONST(CF4_PI,       (float)M_PI) 
DEF_CONST(CF4_0_54,        0.54f) 
DEF_CONST(CF4_0_46,        0.46f) 
DEF_CONST(CF4_2PI,(float)(2.0f*M_PI))
DEF_CONST(CF4_4PI,(float)(4.0f*M_PI))
DEF_CONST(CF4_M1,       -1.0f) 

DEF_CONSTI(CI4_0,        0) 
DEF_CONSTI(CI4_M1,       -1) 

static const int neg_mask = 0xFFFFFFFF;
static const int abs_mask = 0x7FFFFFFF;
static const int sign_mask = 0x80000000;
static const v4sf sign_mask_128 = _mm_set_ps1(((float*)&sign_mask)[0]);
static const v4sf abs_mask_128 = _mm_set_ps1(((float*)&abs_mask)[0]);
static const v4sf neg_mask_128 = _mm_set_ps1(((float*)&neg_mask)[0]);


#ifdef _USE_TRIG_TABLES_  
DEF_CONST(CF4_TTS,(float)trig_table_sz)

float sin_table[trig_table_sz];
float cos_table[trig_table_sz];

void init_trig_tables_v4sf()
{
    float dphi = 4.0f*M_PI/(trig_table_sz-1);
    for(int kk=0;kk<trig_table_sz;kk++)
    {
        cos_table[kk] = std::cos(kk*dphi);
        sin_table[kk] = std::sin(kk*dphi);
    }
}

v4sf v4sf_floor(v4sf x);

v4sf v4sf_fmod(const v4sf x,const v4sf y) 
{ 
    v4sf a;
    a = _mm_div_ps(x,y);
    return _mm_mul_ps(_mm_sub_ps(a,v4sf_floor(a)),y); 
}

#endif


v4sf v4sf_floor(v4sf x)
{
    /*
      movaps       xmm0,   [float_value]
      cvttps2dq    xmm1,   xmm0
      psrld        xmm0,   31
      psubd        xmm1,   xmm0
      cvtsq2ps     xmm0,   xmm1
    */    
    vec4 _xmm0,_xmm1;    
    
    _xmm0.m128f = x;
    _xmm1.m128i = _mm_cvttps_epi32(_xmm0.m128f);
    _xmm0.m128i = _mm_srli_epi32 (_xmm0.m128i,31);    
    _xmm1.m128i = _mm_sub_epi32(_xmm1.m128i,_xmm0.m128i);    
        
    return _mm_cvtepi32_ps(_xmm1.m128i);
}

v4si v4sf_floori(v4sf x)
{
    /*
      movaps       xmm0,   [float_value]
      cvttps2dq    xmm1,   xmm0
      psrld        xmm0,   31
      psubd        xmm1,   xmm0
      cvtsq2ps     xmm0,   xmm1
    */    
    vec4 _xmm0,_xmm1;    
    
    _xmm0.m128f = x;
    _xmm1.m128i = _mm_cvttps_epi32(_xmm0.m128f);
    _xmm0.m128i = _mm_srli_epi32 (_xmm0.m128i,31);    
    _xmm1.m128i = _mm_sub_epi32(_xmm1.m128i,_xmm0.m128i);    
        
    return _xmm1.m128i;
}

v4sf apod_rect_v4sf(v4sf)
{
    return CF4_1;
}
 
v4sf apod_hamming_v4sf(v4sf x)
{    
    return _mm_add_ps(CF4_0_54,_mm_mul_ps(CF4_0_46,cos_ps(_mm_mul_ps(CF4_PI,x))));
}

v4sf xx_func0_v4sf(const v4sf x,const v4sf xd_sz)
{
    return x;
}

v4sf xx_func1_v4sf(const v4sf x,const v4sf xd_sz)
{
    return _mm_sub_ps(x,xd_sz);
}

void planewave_beamforming_rf_v4sf(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_data,const int ch_data_sz[2],float* rf_bf_data,BF_Params params,const float* apodization,const float tx_apodization)
{        
    v4si tix;
    v4sf z,x,x0,rel_ap_sz,ti,ti_frac;
    v4sf x_diff,x_ch;
    
    vec4 apod_mask,tmp,ti_sample,ch_ix,tmp_bf,out_of_rng_mask,apod;
    vec4 chx0,chx1;
    
    int ix_out,tmp32;

    const v4sf ap_lim = _mm_set1_ps(1/params.FN/2.0f);
    
    const v4sf fs_in = _mm_set1_ps(params.fs_in);
    const v4sf fs_out = _mm_set1_ps(params.fs_out);
    const v4sf c0 = _mm_set1_ps(params.c0);
    const v4sf t0 = _mm_set1_ps(params.t0);                

    const v4sf ctxa = cos_ps(_mm_set1_ps(params.tx_angle));
    const v4sf stxa = sin_ps(_mm_set1_ps(params.tx_angle));    
    
    const v4sf crxa = cos_ps(_mm_set1_ps(params.rx_angle));
    const v4sf srxa = sin_ps(_mm_set1_ps(params.rx_angle));       
    
    const v4sf tx_apod = _mm_set1_ps(tx_apodization);
    const int Nch = ch_data_sz[1];
    const int Ns = ch_data_sz[0];
    
    float* out_ptr;
    
#ifdef _HAS_SSE4_    
    const v4si lim = _mm_set1_epi32(Ns-2);
    const v4si lim1 = _mm_set1_epi32(Ns-1);
#else
    const v4sf lim = _mm_set1_ps((float)(Ns-2));    
    const v4sf lim1 = _mm_set1_ps((float)(Ns-1));
#endif            
    
    //const v4sf xd_sz = _mm_set1_ps((Nch-1)*params.dx);
    //XX_FUNCTION xx_func = xx_func0_v4sf;
    //if(params.tx_angle < 0.0f)
    //    xx_func = xx_func1_v4sf;
    
    int offset = 0;
    if(Nr%4 != 0)
        offset = 4;        
    
    int channel_stride = Nx*Nr;
    int apod_ix_offset = 0;
    
    for(int ix=0;ix<Nx;ix++)
    {
        x0 = _mm_set1_ps(x_coords[ix]);        

        for(int iz=0;iz<(Nr-offset);iz+=4)
        {
            ix_out = ix*Nr + iz;
            out_ptr = rf_bf_data + ix_out;
                        
            tmp_bf.m128f = _mm_set1_ps(0.0f);
                                    
            z = _mm_set_ps(r_coords[iz+3],r_coords[iz+2],r_coords[iz+1],r_coords[iz]);
            x = _mm_add_ps(x0,_mm_mul_ps(srxa,z));                                    
            z = _mm_mul_ps(z,crxa);
            
            apod_ix_offset = Nr*ix + iz;
            
            for(int ich=0;ich<Nch;ich++)
            {                    
                x_ch = _mm_set1_ps(ich*params.dx - (((float)(Nch)-1.0f)/2.0f)*params.dx);

                //rel_ap_sz = _mm_div_ps(_mm_and_ps(_mm_sub_ps(x_ch,x0),abs_mask_128),_mm_mul_ps(z,ap_lim));
                //apod_mask.m128f = _mm_cmplt_ps(rel_ap_sz,CF4_1);
                
                // If non of the xd elements are within the aperture, just continue
                int ii = channel_stride*ich + apod_ix_offset;
                apod.m128f = _mm_mul_ps(tx_apod,_mm_set_ps(apodization[ii+3],apodization[ii+1],apodization[ii+2],apodization[ii]));
                
                if(!(apod.f[0] > 0.0f | apod.f[1] > 0.0f | apod.f[2] > 0.0f | apod.f[3] > 0.0f))
                    continue;
                
                x_diff = _mm_sub_ps(x_ch,x);
                
                tmp.m128f = _mm_add_ps(
                        _mm_add_ps(_mm_mul_ps(ctxa,z),_mm_mul_ps(stxa,x)), //_mm_add_ps(_mm_mul_ps(ctxa,z),_mm_mul_ps(stxa,xx_func(x,xd_sz))),
                        _mm_sqrt_ps(_mm_add_ps(_mm_mul_ps(z,z),_mm_mul_ps(x_diff,x_diff)))
                        );
                
                ti = _mm_mul_ps(fs_in,_mm_sub_ps(_mm_div_ps(tmp.m128f,c0),t0));
                
                // ti_sample = std::min(std::max(std::floor((int)ti),0.0f),Ns-2.0f);                                                
                //ti_sample.m128i = v4sf_floori(ti);
#ifdef _HAS_SSE4_
                ti_sample.m128i = _mm_cvttps_epi32(ti);                   
                out_of_rng_mask.m128i = _mm_and_si128(_mm_cmplt_epi32(CI4_M1,ti_sample.m128i),_mm_cmplt_epi32(ti_sample.m128i,lim1));                
                
                ti_sample.m128i = _mm_min_epi32(_mm_max_epi32(ti_sample.m128i,CI4_0),lim);                                
                //ti_frac = ti - ti_sample;
                ti_frac = _mm_sub_ps(ti,_mm_cvtepi32_ps(ti_sample.m128i));                
#else
                ti_sample.m128f = v4sf_floor(ti);
                out_of_rng_mask.m128f = _mm_and_ps(_mm_cmplt_ps(CF4_M1,ti_sample.m128f),_mm_cmplt_ps(ti_sample.m128f,lim1));                
                
                ti_sample.m128f = _mm_min_ps(_mm_max_ps(ti_sample.m128f,CF4_0),lim);                                                
                //ti_frac = ti - ti_sample;
                ti_frac = _mm_sub_ps(ti,ti_sample.m128f);                
                ti_sample.m128i = _mm_cvttps_epi32(ti_sample.m128f);                                
#endif                                
                
                //ch_ix = ich*Ns + tix;
                tmp32 = ich*Ns;
                ch_ix.i[0] = tmp32 + ti_sample.i[0];
                ch_ix.i[1] = tmp32 + ti_sample.i[1];
                ch_ix.i[2] = tmp32 + ti_sample.i[2];
                ch_ix.i[3] = tmp32 + ti_sample.i[3];                                
                
                chx0.f[0] = ch_data[ch_ix.i[0]];
                chx0.f[1] = ch_data[ch_ix.i[1]];
                chx0.f[2] = ch_data[ch_ix.i[2]];
                chx0.f[3] = ch_data[ch_ix.i[3]];
                
                chx1.f[0] = ch_data[ch_ix.i[0]+1];
                chx1.f[1] = ch_data[ch_ix.i[1]+1];
                chx1.f[2] = ch_data[ch_ix.i[2]+1];
                chx1.f[3] = ch_data[ch_ix.i[3]+1];                                
                
//                 tmp.m128f = chx0.m128f;
//                 mexPrintf("%2.5f %2.5f %2.5f %2.5f\n",ch_data[ch_ix.i[0]],tmp.f[1],tmp.f[2],tmp.f[3]);
                //tmpOutRe = ch_data[ch_ix] + (ch_data[ch_ix + 1] - ch_data[ch_ix])*ti_frac;
                tmp.m128f = _mm_add_ps(chx0.m128f,_mm_mul_ps(_mm_sub_ps(chx1.m128f,chx0.m128f),ti_frac));                
                               
                // rf_bf_data[ix_out] += apod*tmpOutRe;                
                //apod_mask.m128f = _mm_and_ps(apod_mask.m128f,out_of_rng_mask.m128f);
                apod_mask.m128f = out_of_rng_mask.m128f;
                tmp_bf.m128f = _mm_add_ps(tmp_bf.m128f,_mm_and_ps(apod_mask.m128f,_mm_mul_ps(apod.m128f,tmp.m128f)));                
            }                  
            
            //_mm_store_ps(out_ptr,tmp_bf);
            out_ptr[0] = tmp_bf.f[0];
            out_ptr[1] = tmp_bf.f[1];
            out_ptr[2] = tmp_bf.f[2];
            out_ptr[3] = tmp_bf.f[3];
        }        
    }          
}

void planewave_beamforming_iq_v4sf(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_dataRe,float* ch_dataIm,const int ch_data_sz[2],float* iq_bf_dataRe,float* iq_bf_dataIm,BF_Params params,const float* apodization,const float tx_apodization)
{
    v4si tix;
    v4sf z,x,x0,rel_ap_sz,ti,ti_frac;
    v4sf x_diff,x_ch,arg;
    
    vec4 apod_mask,tmpRe,tmpIm,ti_sample,ch_ix,arg_ix,out_of_rng_mask,apod;
    vec4 chx0R,chx1R,chx0I,chx1I;
    
    int ix_out,tmp32;

    const v4sf ap_lim = _mm_set1_ps(1/params.FN/2.0f);
        
    const v4sf fs_in = _mm_set1_ps(params.fs_in);
    const v4sf fs_out = _mm_set1_ps(params.fs_out);
    const v4sf c0 = _mm_set1_ps(params.c0);
    const v4sf t0 = _mm_set1_ps(params.t0);
            
    const v4sf w_demod = _mm_set1_ps(2.0f*M_PI*params.f_demod/params.fs_in);    

    const v4sf ctxa = cos_ps(_mm_set1_ps(params.tx_angle));
    const v4sf stxa = sin_ps(_mm_set1_ps(params.tx_angle));    

    vec4 tmp1;
    const v4sf crxa = cos_ps(_mm_set1_ps(params.rx_angle));
    const v4sf srxa = sin_ps(_mm_set1_ps(params.rx_angle));        
    
    const v4sf tx_apod = _mm_set1_ps(tx_apodization);
    
    float dec = params.fs_in/params.fs_out;
    
    v4sf cexp,sexp;            
    
    const int Nch = ch_data_sz[1];
    const int Ns = ch_data_sz[0];
    
    float *out_ptrRe,*out_ptrIm;
    vec4 tmp_bfRe,tmp_bfIm;
#ifdef _HAS_SSE4_    
    const v4si lim = _mm_set1_epi32(Ns-2);
    const v4si lim1 = _mm_set1_epi32(Ns-1);
#else
    const v4sf lim = _mm_set1_ps((float)(Ns-2));    
    const v4sf lim1 = _mm_set1_ps((float)(Ns-1));
#endif      
    
    //const v4sf xd_sz = _mm_set1_ps((Nch-1)*params.dx - (((float)(Nch)-1.0f)/2.0f)*params.dx);
    
    //XX_FUNCTION xx_func = xx_func0_v4sf;
    //if(params.tx_angle < 0.0f)
        //xx_func = xx_func1_v4sf;
    
    int offset = 0;
    if(Nr%4 != 0)
        offset = 4;        
    
    int channel_stride = Nx*Nr;
    int apod_ix_offset = 0;
    
    for(int ix=0;ix<Nx;ix++)
    {        
        x0 = _mm_set1_ps(x_coords[ix]);        
        
        for(int iz=0;iz<(Nr-offset);iz+=4)
        {
            ix_out = ix*Nr + iz;
            
            out_ptrRe = iq_bf_dataRe + ix_out;
            out_ptrIm = iq_bf_dataIm + ix_out;                        
            
            tmp_bfRe.m128f = _mm_set1_ps(0.0f);
            tmp_bfIm.m128f = _mm_set1_ps(0.0f);
            
            z = _mm_set_ps(r_coords[iz+3],r_coords[iz+2],r_coords[iz+1],r_coords[iz]);
            x = _mm_add_ps(x0,_mm_mul_ps(srxa,z));                                    
            z = _mm_mul_ps(z,crxa);
                        
            apod_ix_offset = Nr*ix + iz;
            for(int ich=0;ich<Nch;ich++)
            {                    
                x_ch = _mm_set1_ps(ich*params.dx - (((float)(Nch)-1.0f)/2.0f)*params.dx);
                
                //rel_ap_sz = _mm_div_ps(_mm_and_ps(_mm_sub_ps(x_ch,x0),abs_mask_128),_mm_mul_ps(z,ap_lim));
                //apod_mask.m128f = _mm_cmplt_ps(rel_ap_sz,CF4_1);
                
                int ii = channel_stride*ich + apod_ix_offset;
                apod.m128f = _mm_mul_ps(tx_apod,_mm_set_ps(apodization[ii+3],apodization[ii+1],apodization[ii+2],apodization[ii]));
                // If non of the xd elements are within the aperture, just continue
                if(!(apod.f[0] > 0.0f | apod.f[1] > 0.0f | apod.f[2] > 0.0f | apod.f[3] > 0.0f))
                    continue;                                                
                
                x_diff = _mm_sub_ps(x_ch,x);
                
                tmpRe.m128f = _mm_add_ps(
                        _mm_add_ps(_mm_mul_ps(ctxa,z),_mm_mul_ps(stxa,x)), // _mm_add_ps(_mm_mul_ps(ctxa,z),_mm_mul_ps(stxa,xx_func(x,xd_sz))),
                        _mm_sqrt_ps(_mm_add_ps(_mm_mul_ps(z,z),_mm_mul_ps(x_diff,x_diff)))
                        );
                ti = _mm_mul_ps(fs_in,_mm_sub_ps(_mm_div_ps(tmpRe.m128f,c0),t0));
                
                // ti_sample = std::min(std::max(std::floor((int)ti),0.0f),Ns-2.0f);                                                
                //ti_sample.m128i = v4sf_floori(ti);
#ifdef _HAS_SSE4_
                ti_sample.m128i = _mm_cvttps_epi32(ti);                 
                out_of_rng_mask.m128i = _mm_and_si128(_mm_cmplt_epi32(CI4_M1,ti_sample.m128i),_mm_cmplt_epi32(ti_sample.m128i,lim1));                        
                ti_sample.m128i = _mm_min_epi32(_mm_max_epi32(ti_sample.m128i,CI4_0),lim);                                
                //ti_frac = ti - ti_sample;
                ti_frac = _mm_sub_ps(ti,_mm_cvtepi32_ps(ti_sample.m128i));                
#else
                ti_sample.m128f = v4sf_floor(ti);
                out_of_rng_mask.m128f = _mm_and_ps(_mm_cmplt_ps(CF4_M1,ti_sample.m128f),_mm_cmplt_ps(ti_sample.m128f,lim1));                
                ti_sample.m128f = _mm_min_ps(_mm_max_ps(ti_sample.m128f,CF4_0),lim);                                                
                //ti_frac = ti - ti_sample;
                ti_frac = _mm_sub_ps(ti,ti_sample.m128f);                
                ti_sample.m128i = _mm_cvttps_epi32(ti_sample.m128f);                                
#endif                
                //ch_ix = ich*Ns + tix;
                tmp32 = ich*Ns;
                ch_ix.i[0] = tmp32 + ti_sample.i[0];
                ch_ix.i[1] = tmp32 + ti_sample.i[1];
                ch_ix.i[2] = tmp32 + ti_sample.i[2];
                ch_ix.i[3] = tmp32 + ti_sample.i[3];                                
                
                chx0R.f[0] = ch_dataRe[ch_ix.i[0]];
                chx0R.f[1] = ch_dataRe[ch_ix.i[1]];
                chx0R.f[2] = ch_dataRe[ch_ix.i[2]];
                chx0R.f[3] = ch_dataRe[ch_ix.i[3]];
                
                chx1R.f[0] = ch_dataRe[ch_ix.i[0]+1];
                chx1R.f[1] = ch_dataRe[ch_ix.i[1]+1];
                chx1R.f[2] = ch_dataRe[ch_ix.i[2]+1];
                chx1R.f[3] = ch_dataRe[ch_ix.i[3]+1];                                
                
                chx0I.f[0] = ch_dataIm[ch_ix.i[0]];
                chx0I.f[1] = ch_dataIm[ch_ix.i[1]];
                chx0I.f[2] = ch_dataIm[ch_ix.i[2]];
                chx0I.f[3] = ch_dataIm[ch_ix.i[3]];
                
                chx1I.f[0] = ch_dataIm[ch_ix.i[0]+1];
                chx1I.f[1] = ch_dataIm[ch_ix.i[1]+1];
                chx1I.f[2] = ch_dataIm[ch_ix.i[2]+1];
                chx1I.f[3] = ch_dataIm[ch_ix.i[3]+1];                                
                
                //tmpOutRe = ch_data[ch_ix] + (ch_data[ch_ix + 1] - ch_data[ch_ix])*ti_frac;
                tmpRe.m128f = _mm_add_ps(chx0R.m128f,_mm_mul_ps(_mm_sub_ps(chx1R.m128f,chx0R.m128f),ti_frac));
                tmpIm.m128f = _mm_add_ps(chx0I.m128f,_mm_mul_ps(_mm_sub_ps(chx1I.m128f,chx0I.m128f),ti_frac));
                
                arg = _mm_mul_ps(w_demod,_mm_sub_ps(ti,_mm_set_ps(dec*(iz+3),dec*(iz+2),dec*(iz+1),dec*(iz+0))));
 #ifdef _USE_TRIG_TABLES_
                arg = _mm_add_ps(v4sf_fmod(arg,CF4_2PI),CF4_2PI);
                arg_ix.m128i = _mm_cvttps_epi32(_mm_div_ps(_mm_mul_ps(CF4_TTS,arg),CF4_4PI));
                
                cexp = _mm_set_ps(cos_table[arg_ix.i[3]],
                        cos_table[arg_ix.i[2]],cos_table[arg_ix.i[1]],cos_table[arg_ix.i[0]]);
                sexp = _mm_set_ps(sin_table[arg_ix.i[3]],
                        sin_table[arg_ix.i[2]],sin_table[arg_ix.i[1]],sin_table[arg_ix.i[0]]);
#else
                cexp = cos_ps(arg);
                sexp = sin_ps(arg);
#endif
//                 
//                 tmpRe.m128f = arg;
//                 mexPrintf("%2.5f %2.5f %2.5f %2.5f\n",tmpRe.f[0],tmpRe.f[1],tmpRe.f[2],tmpRe.f[3]);
                                                    
                //apod_mask.m128f = _mm_and_ps(apod_mask.m128f,out_of_rng_mask.m128f);
                apod_mask.m128f = out_of_rng_mask.m128f;
                
                tmp_bfIm.m128f = _mm_add_ps(tmp_bfIm.m128f,
                        _mm_and_ps(apod_mask.m128f,_mm_mul_ps(apod.m128f,_mm_add_ps(_mm_mul_ps(tmpIm.m128f,cexp),_mm_mul_ps(tmpRe.m128f,sexp))))
                        );                
                
                tmpRe.m128f = _mm_sub_ps(_mm_mul_ps(tmpRe.m128f,cexp),_mm_mul_ps(tmpIm.m128f,sexp));
                tmp_bfRe.m128f = _mm_add_ps(tmp_bfRe.m128f,_mm_and_ps(apod_mask.m128f,_mm_mul_ps(apod.m128f,tmpRe.m128f)));                                
            }                                  
            out_ptrRe[0] = tmp_bfRe.f[0];
            out_ptrRe[1] = tmp_bfRe.f[1];
            out_ptrRe[2] = tmp_bfRe.f[2];
            out_ptrRe[3] = tmp_bfRe.f[3];
            
            out_ptrIm[0] = tmp_bfIm.f[0];
            out_ptrIm[1] = tmp_bfIm.f[1];
            out_ptrIm[2] = tmp_bfIm.f[2];
            out_ptrIm[3] = tmp_bfIm.f[3];
            
            //_mm_store_ps(out_ptrRe,tmp_bfRe);
            //_mm_store_ps(out_ptrIm,tmp_bfIm);
        }        
    }
}
