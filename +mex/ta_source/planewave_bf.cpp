    #include "planewave_bf.hpp"
#include <cmath>
#include <algorithm>
#include <mex.h>

#include <iostream>
#include <fstream>
    
const float M_2PIf = 2.0f*M_PI;
const float M_4PIf = 4.0f*M_PI;

typedef float (*APOD_FUNCTION)(float);
typedef float (*XX_FUNCTION)(float,float);
#ifdef _USE_TRIG_TABLES_  
float sin_table[trig_table_sz];
float cos_table[trig_table_sz];

void init_trig_tables()
{
    float dphi = M_4PIf/(trig_table_sz-1);
    for(int kk=0;kk<trig_table_sz;kk++)
    {
        cos_table[kk] = std::cos(kk*dphi);
        sin_table[kk] = std::sin(kk*dphi);
    }
}

float mfmod(const float x,const float y) 
{ 
    float a; return ((a=x/y)-(int)a)*y; 
}

#endif

float apod_rect(float)
{
    return 1.0f;
}

float apod_hamming(float x)
{
    return 0.54f + 0.46f*std::cos(M_PI*x);
}

float xx_func0(const float x,const float xd_sz)
{
    return x;
}

float xx_func1(const float x,const float xd_sz)
{
    return x - xd_sz;
}

float limit(const float x,const float lo,const float hi,float &amp)
{        
    if(x > hi)
    {
        amp = 0.0f;
        return hi;
    }
    else if(x < lo)
    {
        amp = 0.0f;
        return lo;
    }
    
    return x;
}

void planewave_beamforming_rf(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_data,const int ch_data_sz[2],float* rf_bf_data,BF_Params params,const float* apod_mask)
{
    float ti,ti_sample,ti_frac,z_lim,tmpOutRe;    
    float z,x,x0,x_diff,x_ch;

    int tix, ix_out,ch_ix;

    const float ap_lim = 1/params.FN/2.0f;
    float apod = 1.0f;
    const float w_demod = 2.0f*M_PI*params.f_demod/params.fs_in;    

    const float ctxa = std::cos(params.tx_angle);
    const float stxa = std::sin(params.tx_angle);

    const float crxa = std::cos(params.rx_angle);
    const float srxa = std::sin(params.rx_angle);              
    
    int Nch = ch_data_sz[1];
    int Ns = ch_data_sz[0];        
    
    const float xd_sz = (Nch-1)*params.dx;
    XX_FUNCTION xx_func = xx_func0;
    if(params.tx_angle < 0.0f)
        xx_func = xx_func1;       
    
    const int channel_stride = Nx*Nr;
    int apod_ix_offset = 0;            
    for(int ix=0;ix<Nx;ix++)
    {        
        x0 = x_coords[ix];        
        for(int iz=0;iz<Nr;iz++)
        {
            ix_out = ix*Nr + iz;
            rf_bf_data[ix_out] = 0;

            z = crxa*r_coords[iz];
            x = x0 + srxa*r_coords[iz];                                    
            apod_ix_offset = Nr*ix + (iz-0);
            for(int ich=0;ich<Nch;ich++)
            {                                                
                x_ch = ich*params.dx;
                
                //rel_ap_sz = std::abs(x_ch - x0)/z/ap_lim;
                apod = tx_apodization*apod_mask[channel_stride*ich + apod_ix_offset];
                
                if(apod == 0.0f)
                    continue;                
                
                x_diff = x_ch - x;
                ti = params.fs_in*((ctxa*z + stxa*xx_func(x,xd_sz) + std::sqrt(z*z + x_diff*x_diff))/params.c0 - params.t0);                                
                
                //ti_sample = std::min(std::max(std::floor(ti),0.0f),Ns-2.0f);                
                ti_sample = limit(std::floor(ti),0.0f,Ns-2.0f,apod);
                
                ti_frac = ti - ti_sample;
                tix = (int)ti_sample;
                                
                ch_ix = (ich*Ns + tix);
                tmpOutRe = ch_data[ch_ix] + (ch_data[ch_ix + 1] - ch_data[ch_ix])*ti_frac;
     
                rf_bf_data[ix_out] += apod*tmpOutRe;                
            }            
        }
    }      
}

void planewave_beamforming_iq(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_dataRe,float* ch_dataIm,const int ch_data_sz[2],float* iq_bf_dataRe,float* iq_bf_dataIm,BF_Params params,const float* apod_mask)
{
    float ti,ti_sample,ti_frac,z_lim,tmpOutRe,tmpOutIm;    
    float z,x,x0,x_diff,x_ch,arg;

    int tix, ix_out,ch_ix,arg_indx;

    float ap_lim = 1/params.FN/2.0f;
    float apod = 1.0f;
    const float w_demod = 2.0f*M_PI*params.f_demod/params.fs_in;    

    float ctxa = std::cos(params.tx_angle);
    float stxa = std::sin(params.tx_angle);

    float crxa = std::cos(params.rx_angle);
    float srxa = std::sin(params.rx_angle);      
    
    float cexp,sexp;            
    
    int Nch = ch_data_sz[1];
    int Ns = ch_data_sz[0];        
    
    const float xd_sz = (Nch-1)*params.dx;
    XX_FUNCTION xx_func = xx_func0;
    if(params.tx_angle < 0.0f)
        xx_func = xx_func1;       
    
    float dec = params.fs_in/params.fs_out;
    
    int channel_stride = Nx*Nr;
    int apod_ix_offset = 0;
    for(int ix=0;ix<Nx;ix++)
    {                
        x0 = x_coords[ix];        
        
        for(int iz=0;iz<Nr;iz++)
        {
            ix_out = ix*Nr + iz;
            
            iq_bf_dataRe[ix_out] = 0.0f;
            iq_bf_dataIm[ix_out] = 0.0f;
            
            z = crxa*r_coords[iz];
            x = x0 + srxa*r_coords[iz];                                    

            apod_ix_offset = Nr*ix + iz;
            for(int ich=0;ich<Nch;ich++)
            {                                                
                x_ch = ich*params.dx;                

                apod = tx_apodization*apod_mask[channel_stride*ich + apod_ix_offset];
                
                if(apod == 0.0f)
                    continue;                
                
                x_diff = x_ch - x;
                ti = params.fs_in*((ctxa*z + stxa*xx_func(x,xd_sz) + std::sqrt(z*z + x_diff*x_diff))/params.c0 - params.t0);                                
                                
                //ti_sample = std::min(std::max(std::floor(ti),0.0f),Ns-2.0f);                
                ti_sample = limit(std::floor(ti),0.0f,Ns-2.0f,apod);
                
                ti_frac = ti - ti_sample;
                tix = (int)ti_sample;
                ch_ix = ich*Ns + tix;
                
                tmpOutRe = ch_dataRe[ch_ix] + (ch_dataRe[ch_ix+1] - ch_dataRe[ch_ix])*ti_frac;
                tmpOutIm = ch_dataIm[ch_ix] + (ch_dataIm[ch_ix+1] - ch_dataIm[ch_ix])*ti_frac;   
                                                
                arg = w_demod*(ti - dec*iz);
                
#ifdef _USE_TRIG_TABLES_
                arg = mfmod(arg,M_2PIf) + M_2PIf;
                arg_indx = (int)(trig_table_sz*arg/M_4PIf);
                cexp = cos_table[arg_indx];
                sexp = sin_table[arg_indx];
#else
                cexp = std::cos(arg);
                sexp = std::sin(arg);
#endif
                iq_bf_dataIm[ix_out] += apod*(tmpOutIm*cexp + tmpOutRe*sexp);                    
                tmpOutRe = tmpOutRe*cexp - tmpOutIm*sexp;                    
                        
                iq_bf_dataRe[ix_out] += apod*tmpOutRe;                
            }            
        }
    }
}