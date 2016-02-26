#include <mex.h>
#include <xmmintrin.h>
#include <cmath>
#include "common.hpp"
#include "planewave_bf.hpp"

#if defined(_WIN_) && defined (_PAR_FOR_) 
    #include <ppl.h> // requires VS2010+
#endif
        
#if defined (_UNIX_) && defined (_PAR_FOR_) 
    #include <omp.h>        
#endif        
        
static float BF_PI = 3.141592653589793f;
static bool trig_table_initialized = false;

template<class T>
        T my_min(const T& a,const T& b)
{
    return (a < b ? a : b);
}

template<class T>
        T my_max(const T& a,const T& b)
{
    return (a < b ? b : a);
}
extern"C" {
    
    void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
      const mxArray *prhs[])
    {
#ifdef _USE_TRIG_TABLES_        
        if(!trig_table_initialized)
        {
            mexPrintf("Initializing trig tables\n");
            
#ifdef _USE_V4SF_            
            init_trig_tables_v4sf();
#else
            init_trig_tables();
#endif
            trig_table_initialized = true;
        }
#endif        
        if(nrhs < 5)
        {
            mexErrMsgTxt("Not enough input arguments");
        }

        if(!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[2]) || !mxIsSingle(prhs[3]) || !mxIsSingle(prhs[4]))
        {
            mexErrMsgTxt("Input data must be of type single");
        }
/*
        if(nrhs > 5)
        {
            if(!mxIsSingle(prhs[5]))
            {
                mexErrMsgTxt("Transmit apodization must be of type single");
            }  
        }
*/
        const mwSize* iDims = mxGetDimensions(prhs[0]);
        int iNumberOfDimensions = mxGetNumberOfDimensions(prhs[0]);
        if(iNumberOfDimensions < 2)
            mexErrMsgTxt("Channel data must have more than 2 dimensions");
        
        const int Ns = iDims[0];
        const int Nch = iDims[1];
        const int Nf = (iNumberOfDimensions > 2) ? iDims[2] : 1;

        float* iXcoords = (float*)(mxGetPr(prhs[2]));
        float* iZcoords = (float*)(mxGetPr(prhs[3]));

        int Nx = mxGetN(prhs[2])*mxGetM(prhs[2]);
        int Nz = mxGetN(prhs[3])*mxGetM(prhs[3]);

        if(nlhs == 0)
            return;

        mwSize iOutDims[3] = {Nz,Nx,Nf};

        float *iOutDataRe = NULL,*iOutDataIm = NULL;    
        float *iChDataRe = NULL,*iChDataIm = NULL;
        float *iApodization = NULL;

        bool isComplex = mxIsComplex(prhs[0]);
        
        if(isComplex)
        {
            plhs[0] = mxCreateNumericArray(3,iOutDims,mxSINGLE_CLASS,mxCOMPLEX);

            iOutDataRe = (float*)(mxGetPr(plhs[0]));
            iOutDataIm = (float*)(mxGetPi(plhs[0]));

            iChDataRe = (float*)(mxGetPr(prhs[0]));
            iChDataIm = (float*)(mxGetPi(prhs[0]));
        }
        else
        {
            plhs[0] = mxCreateNumericArray(3,iOutDims,mxSINGLE_CLASS,mxREAL);
            iOutDataRe = (float*)(mxGetPr(plhs[0]));        
            iChDataRe = (float*)(mxGetPr(prhs[0]));        
        }                                      
        
        float fs_in,fs_out,t0,c0,dx,FN,f_demod = 0.0f,rx_angle;
        const mxArray* iParams = prhs[1];
        if(!GetScalar(mxGetField(iParams,0,"fs_in"),fs_in))
        {
            mexErrMsgTxt("fs_in must be a field in params struct.");
        }

        if(!GetScalar(mxGetField(iParams,0,"fs_out"),fs_out))
        {
            mexErrMsgTxt("fs_out must be a field in params struct.");
        }
        
        if(!GetScalar(mxGetField(iParams,0,"t0"),t0))
        {
            mexErrMsgTxt("t0 must be a field in params struct.");
        }       

        if(!GetScalar(mxGetField(iParams,0,"c"),c0))
        {
            mexErrMsgTxt("c0 must be a field in params struct.");
        }

        if(!GetScalar(mxGetField(iParams,0,"pitch"),dx))
        {
            mexErrMsgTxt("pitch must be a field in params struct.");
        }

        if(!GetScalar(mxGetField(iParams,0,"rx_angle"),rx_angle))
        {
            mexErrMsgTxt("dx must be a field in params struct.");
        }
        
        if(!GetScalar(mxGetField(iParams,0,"FN"),FN))
        {
            mexErrMsgTxt("dx must be a field in params struct.");
        }

        if(isComplex)
        {
            if(!GetScalar(mxGetField(iParams,0,"f_demod"),f_demod))
            {
                mexErrMsgTxt("f_demod must be a field in params struct.");
            }
        }
        
        if(mxGetField(iParams,0,"tx_angle") == NULL)
        {
            mexErrMsgTxt("tx_angle must be a field in params struct.");
        }
        
        mxArray* tx_angles = mxGetField(iParams,0,"tx_angle");
        const int num_tx_angles = mxGetN(tx_angles)*mxGetM(tx_angles);
        if(num_tx_angles > 1 && num_tx_angles < Nf)
        {
            mexErrMsgTxt("Number of tx angles must be either 1 or at least size(ch_data,3).");
        }
        
        float* iTxApodization = new float[num_tx_angles];
       
        if(nrhs > 5)
        {
            const mxArray* tx_apodizationMx = prhs[5];
            if(mxGetN(prhs[5])*mxGetM(prhs[5]) != num_tx_angles)
            {
                mexErrMsgTxt("Tx apodization must have length(tx_angles) elements");
            }
            
            for(int kk=0;kk<num_tx_angles;kk++)
                iTxApodization[kk] = ((float*)mxGetPr(prhs[5]))[kk];
        }
        else
        {
            for(int kk=0;kk<num_tx_angles;kk++)
                iTxApodization[kk] = 1.0f;
        }
                
#ifndef _PAR_FOR_        
        BF_Params params;
        params.tx_angle = 0;
        params.rx_angle = rx_angle;
        params.dx = dx;
        params.fs_in = fs_in;
        params.fs_out = fs_out;
        params.t0 = t0;
        params.c0 = c0;
        params.FN = FN;
        params.f_demod = f_demod; 
        
        iApodization = (float*)(mxGetPr(prhs[4]));
        const mwSize* iApodMaskSize = mxGetDimensions(prhs[4]);
        iNumberOfDimensions = mxGetNumberOfDimensions(prhs[4]);
        if(iNumberOfDimensions != 3)
        {
            mexErrMsgTxt("Apodization mask must have 3 dimensions");
        }
        
        if(iApodMaskSize[0] != Nz || iApodMaskSize[1] != Nx ||  iApodMaskSize[2] != Nch)
        {
            mexErrMsgTxt("Apodization mask must have be of size [length(z) length(x) size(ch_data,2)].");
        }        
        
#endif        
        const int chDataSz[2] = {Ns,Nch};
        
#if defined (_WIN_) && defined (_PAR_FOR_) 
        //for(int frame_ix = 0;frame_ix < Nf;frame_ix++)
        Concurrency::parallel_for(0,Nf,[&](int frame_ix)
#else
    #if defined (_UNIX_) && defined (_PAR_FOR_) 
        int frame_ix;    
        #pragma omp parallel num_threads(omp_get_num_procs())
        {
        #pragma omp for
        for(frame_ix = 0;frame_ix < Nf;++frame_ix)
            
    #else          
        for(int frame_ix = 0;frame_ix < Nf;frame_ix++)
    #endif
#endif            
        {            
        
#ifdef _PAR_FOR_        
            BF_Params params;
            params.tx_angle = 0;
            params.rx_angle = rx_angle;
            params.dx = dx;
            params.fs_in = fs_in;
            params.fs_out = fs_out;
            params.t0 = t0;
            params.c0 = c0;
            params.FN = FN;
            params.f_demod = f_demod;            

            iApodization = (float*)(mxGetPr(prhs[4]));
            const mwSize* iApodMaskSize = mxGetDimensions(prhs[4]);
            iNumberOfDimensions = mxGetNumberOfDimensions(prhs[4]);
            if(iNumberOfDimensions != 3)
            {
                mexErrMsgTxt("Apodization mask must have 3 dimensions");
            }

            if(iApodMaskSize[0] != Nz || iApodMaskSize[1] != Nx ||  iApodMaskSize[2] != Nch)
            {
                mexErrMsgTxt("Apodization mask must have be of size [length(z) length(x) size(ch_data,3)].");
            }        
#endif                
            switch(mxGetClassID(tx_angles))
            {
                case mxSINGLE_CLASS:                    
                    params.tx_angle = ((float*)mxGetPr(tx_angles))[(num_tx_angles == 1) ? 0 : frame_ix];
                    break;
                case mxDOUBLE_CLASS:
                    params.tx_angle = (float)(mxGetPr(tx_angles)[(num_tx_angles == 1) ? 0 : frame_ix]);
                    break;
                default:
                    mexErrMsgTxt("p.tx_angle must be either single or double");
            }                                    
            float tx_apodization = iTxApodization[(num_tx_angles == 1) ? 0 : frame_ix];
            int in_offset = frame_ix*Ns*Nch;
            int out_offset = frame_ix*Nx*Nz;
            
#ifdef _USE_V4SF_
            if(!isComplex)
            {            
                planewave_beamforming_rf_v4sf(iXcoords,Nx,iZcoords,Nz,iChDataRe+in_offset,chDataSz,iOutDataRe+out_offset,params,iApodization,tx_apodization);            
            }
            else
            {
                planewave_beamforming_iq_v4sf(iXcoords,Nx,iZcoords,Nz,iChDataRe+in_offset,iChDataIm+in_offset,chDataSz,iOutDataRe+out_offset,iOutDataIm+out_offset,params,iApodization,tx_apodization);            
            }
#else
            if(!isComplex)
            {            
                planewave_beamforming_rf(iXcoords,Nx,iZcoords,Nz,iChDataRe+in_offset,chDataSz,iOutDataRe+out_offset,params,iApodization,tx_apodization);            
            }
            else
            {
                planewave_beamforming_iq(iXcoords,Nx,iZcoords,Nz,iChDataRe+in_offset,iChDataIm+in_offset,chDataSz,iOutDataRe+out_offset,iOutDataIm+out_offset,params,iApodization,tx_apodization);            
            }
#endif

#if defined (_WIN_) && defined (_PAR_FOR_)
        });
#else
        }
#if defined (_UNIX_) && defined (_PAR_FOR_)        
    }
#endif
#endif
        delete[] iTxApodization;
        const mwSize iOutDims2[3] = {Nz,Nx};                                        
        if(nlhs > 1)
        {         
            plhs[1] = mxCreateNumericArray(2,iOutDims2,mxSINGLE_CLASS,mxREAL);
            
            float* Xgrid = (float*)mxGetPr(plhs[1]);            
            for(int ix=0;ix<Nx;ix++)
            {
                for(int iz=0;iz<Nz;iz++)
                {
                    Xgrid[ix*Nz + iz] = iXcoords[ix] + std::sin(rx_angle)*iZcoords[iz];
                }
            }
        }        

        if(nlhs > 2)
        {         
            plhs[2] = mxCreateNumericArray(2,iOutDims2,mxSINGLE_CLASS,mxREAL);
            
            float* Zgrid = (float*)mxGetPr(plhs[2]);            
            for(int ix=0;ix<Nx;ix++)
            {
                for(int iz=0;iz<Nz;iz++)
                {
                    Zgrid[ix*Nz + iz] = std::cos(rx_angle)*iZcoords[iz];
                }
            }
        }        
    }      
}
