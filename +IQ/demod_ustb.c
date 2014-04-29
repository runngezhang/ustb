#include "mex.h"
#include <string.h>
#include <ipp.h>

void demodulate_hilbert(const double *rf,const int length,const mwSize cols,const double f_demod,const double fs,mxArray* iq)
{
    int k;
    IppsHilbertSpec_32f32fc* hilbertSpec;
    
    Ipp32f* tmpRf = ippsMalloc_32f(length);
    
    Ipp32f* tmpReal = ippsMalloc_32f(length);
    Ipp32f* tmpImag = ippsMalloc_32f(length);

    Ipp32fc* mix = ippsMalloc_32fc(length);
    
    Ipp32fc* tmpIq = ippsMalloc_32fc(length);
    
    Ipp64f* real = mxGetPr(iq);
    Ipp64f* imag = mxGetPi(iq);
    
    const double *ptrRf = rf,Ts = 1/fs;
    double *ptrReal = real,*ptrImag = imag;
    
    for(k=0;k<length;k++)
    {
        mix[k].re = (Ipp32f)cos(-IPP_2PI*f_demod*k*Ts);
        mix[k].im = (Ipp32f)sin(-IPP_2PI*f_demod*k*Ts);
    }
    
    ippsHilbertInitAlloc_32f32fc(&hilbertSpec,length,ippAlgHintFast);
    
    for(k=0;k<cols;k++)
    {
        ippsConvert_64f32f(ptrRf,tmpRf,length);
        
        ippsHilbert_32f32fc(tmpRf,tmpIq,hilbertSpec);
        
        ippsCplxToReal_32fc(tmpIq,tmpReal,tmpImag,length);
        ippsMul_32fc_I(mix,tmpIq,length);
        
        ippsConvert_32f64f(tmpReal,ptrReal,length);
        ippsConvert_32f64f(tmpImag,ptrImag,length);
        
        ptrRf += length;
        ptrReal += length;
        ptrImag += length;
    }    
    
    /* Free data */
    ippsHilbertFree_32f32fc(hilbertSpec);
    ippsFree(tmpRf);        
    ippsFree(tmpIq);
    ippsFree(tmpImag);
    ippsFree(tmpReal);
    ippsFree(mix);
}

void demodulate_narrowband(const double *rf,const int length,const mwSize cols,const double f_demod,const double fs,mxArray* iq,const double* firCoefs,const int numCoefs)
{     
    int k;
    Ipp64f* cosTable = ippsMalloc_64f(length);
    Ipp64f* sinTable = ippsMalloc_64f(length);
    Ipp64f* tmp1 = ippsMalloc_64f(length);
    Ipp64f* tmp2 = ippsMalloc_64f(length);
    
    Ipp64f* real = mxGetPr(iq);
    Ipp64f* imag = mxGetPi(iq);
    Ipp64f* delayLine = ippsMalloc_64f(numCoefs);    
    
    const double *ptrRf = rf,Ts = 1/fs;
    double *ptrReal = real,*ptrImag = imag;
    IppsFIRState_64f* firState;        
    
    ippsSet_64f(0.0,delayLine,numCoefs);    
    ippsFIRInitAlloc_64f(&firState,firCoefs,numCoefs,delayLine);        
    
    for(k=0;k<length;k++)
    {
        cosTable[k] = cos(-IPP_2PI*f_demod*k*Ts);
        sinTable[k] = sin(-IPP_2PI*f_demod*k*Ts);
    }    
    
    for(k=0;k<cols;k++)
    {
        ippsMul_64f(ptrRf,cosTable,tmp1,length);

        ippsFIR_64f(tmp1,tmp2,length,firState);        
        ippsFIRSetDlyLine_64f(firState,delayLine);
        
        ippsFlip_64f(tmp2,tmp1,length);
        
        ippsFIR_64f(tmp1,tmp2,length,firState);        
        ippsFIRSetDlyLine_64f(firState,delayLine);
        
        ippsFlip_64f(tmp2,ptrReal,length);
        
        ippsMul_64f(ptrRf,sinTable,tmp1,length);
        
        ippsFIR_64f(tmp1,tmp2,length,firState);        
        ippsFIRSetDlyLine_64f(firState,delayLine);
        
        ippsFlip_64f(tmp2,tmp1,length);
        
        ippsFIR_64f(tmp1,tmp2,length,firState);        
        ippsFIRSetDlyLine_64f(firState,delayLine);
        
        ippsFlip_64f(tmp2,ptrImag,length);
        
        ptrRf += length;
        ptrReal += length;
        ptrImag += length;
    }    
        
    ippsFIRFree_64f(firState);
    ippsFree(cosTable);
    ippsFree(sinTable);   
    ippsFree(delayLine);
    ippsFree(tmp1);
    ippsFree(tmp2);    
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    char* method;
    double fs,f_demod;
    double* firCoefs;
    double cutOff;
    
    int length,k,numCoefs,ippFir = 0;
    mwSize ndim,cols;
    mwSize* dims;
    
    if(nrhs < 4)
    {
        mexErrMsgTxt("Too few input arguments");
        return;
    }
    
    fs = mxGetPr(prhs[1])[0];
    f_demod = mxGetPr(prhs[2])[0];
    
    method = mxMalloc(mxGetN(prhs[3])*mxGetM(prhs[3])+1);
    mxGetString(prhs[3],method,20);
    
    length = mxGetM(prhs[0]);

    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    
    if(ndim > 2)
    {
        mexErrMsgTxt("Only works for matrices. Too many dimensions");
        return;
    }
    /*
    mexPrintf("ndim %d\n",ndim);
    mexPrintf("N %d\n",dims[0]);
    mexPrintf("M %d\n",dims[1]);
    mexPrintf("length %d\n",length);
    */
    //plhs[0] = mxCreateNumericMatrix(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);        
    plhs[0] = mxCreateDoubleMatrix(dims[0],dims[1],mxCOMPLEX);    
    
    if(nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        mxGetPr(plhs[1])[0] = 1.0;
    }
    
    cols = 1;
    for(k=1;k<ndim;k++)
    {
        cols *= dims[k];
    }
    
    if(strcmp(method,"hilbert") == 0)
    {        
        demodulate_hilbert(mxGetPr(prhs[0]),length,cols,f_demod,fs,plhs[0]);    
        if(nlhs > 2)
        {
            plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(plhs[2])[0] = 1.0;
        }
    }
    else if(strcmp(method,"narrowband") == 0)
    {
        if(nrhs < 6)
        {
            mexErrMsgTxt("Too few arguments.\n");
            return;
        }
        
        if(mxGetN(prhs[4]) > 1 || mxGetM(prhs[4]) > 1)
        {
            firCoefs = mxGetPr(prhs[4]);
            numCoefs = max(mxGetN(prhs[4]),mxGetM(prhs[4]));            
        }
        else
        {
            numCoefs = mxGetPr(prhs[4])[0];
            cutOff = mxGetPr(prhs[5])[0];
            cutOff = 1e6*cutOff/fs;
            firCoefs = ippsMalloc_64f(numCoefs);            
            ippsFIRGenLowpass_64f(cutOff,firCoefs,numCoefs,ippWinHamming,ippFalse);  
            ippFir = 1;
        }
        
        demodulate_narrowband(mxGetPr(prhs[0]),length,cols,f_demod,fs,plhs[0],firCoefs,numCoefs);


        if(nlhs > 2)
        {
            plhs[2] = mxCreateDoubleMatrix(1,numCoefs,mxREAL);
            ippsCopy_64f(firCoefs,mxGetPr(plhs[2]),numCoefs);
        }
        
        if(ippFir)
        {
            ippsFree(firCoefs);
        }
    }

    mxFree(method);

}