#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>
#include <concurrent_vector.h>
#include <tuple>
#include <valarray>
#include <set>
#if defined(_WIN_) 
    #include <ppl.h>           // Requires VS2010+
#elif defined (_UNIX_) 
    #include <parallel_for.h>  // Requires Intel tbb
#endif    

// compulsory input
#define	M_P         prhs[0]	// channel_data [time, channel, frame]
#define	M_FS		prhs[1] // sampling frequency (Hz)
#define M_T0		prhs[2]	// initial time (s)
#define	M_DELAY	    prhs[3]	// delays [pixel, channel]
#define	M_APO	    prhs[4]	// apodization % phase correction [pixel, channel]
#define	M_FD		prhs[5] // modulation frequency (Hz)

// optional input
#define	M_VERBOSE	prhs[6] // verbose flag

// output
#define	M_D			plhs[0] // delayed data [pixel, channel, frame]

// constants
#define EPS 1e-6
#define PI 3.14159265359

// types
typedef std::vector<float*> vec_p_float; // change to single

//-----------------------------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	///////////////////////////////
	// CHECKING ARGUMENTS
	///////////////////////////////////////
	// number of inputs
	if (nrhs<6) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nrhs", "Too few input arguments");
	if (nlhs>1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nlhs", "Too many output arguments");

	///////////////////////////////////////
	// VERBOSE
	///////////////////////////////////////
	bool verbose = false;
	if (nrhs == 7)  verbose = (*((float*)mxGetData(M_VERBOSE))) > EPS;

	if (verbose) {
		mexPrintf("---------------------------------------------------------------\n");
		mexPrintf(" USTB mex delay-and-sum for multiple waves\n");
		mexPrintf("---------------------------------------------------------------\n");
		mexPrintf(" Single precision\n");
		mexPrintf(" Vers:  1.0.4\n");
		mexPrintf(" Auth:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>\n");
		mexPrintf(" Date:  2017/05/02\n");
		mexPrintf("---------------------------------------------------------------\n");
	}

	// Channel data
	int ndim = (int)mxGetNumberOfDimensions(M_P);
	if (ndim<2 || ndim>4) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown channel data format. Expected from 2 to 4 dimensions: [time, channel, wave, frame]");
	const mwSize* p_dim = mxGetDimensions(M_P);
	int L = (int)p_dim[0];	// number of time samples
	int N = (int)p_dim[1];	// number of channels 
	int W = (int)1; if (ndim > 2) W = (int)p_dim[2];	// number of waves
    int F = (int)1;	if (ndim > 3) F = (int)p_dim[3];	// number of frames
	if (mxIsDouble(M_P)) mexErrMsgTxt("The channel data should be single precision");    
	
    if (verbose) {
		mexPrintf("Time Samples                    %i\n", L);
		mexPrintf("Channels						%i\n", N);
		mexPrintf("Waves                            %i\n", W);        
		mexPrintf("Frames							%i\n", F);
	}

    // Delays
	ndim = (int)mxGetNumberOfDimensions(M_DELAY);
	if (ndim<2 || ndim>3) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown delay format. Expected 2 or 3 dimensions: [pixel, channel, waves]");
    p_dim = mxGetDimensions(M_DELAY);
	int P = (int)p_dim[0];	// number of pixels
	int N_check = (int)p_dim[1];	// number of channels 
    int W_check = (int)1; if (ndim == 3) W_check = (int)p_dim[2];	// number of waves
    if (N!=N_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of channels in channel data & requested delays do not match.");
	if (W!=W_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of waves in channel data & requested delays do not match.");
    if (mxIsDouble(M_DELAY)) mexErrMsgTxt("The delay should be single precision");    
	
	if (verbose) {
		mexPrintf("Pixels							%i\n", P);
	}

    // Apodization
	ndim = (int)mxGetNumberOfDimensions(M_APO);
	if (ndim<2 || ndim>3) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown apodization format. Expected 2 o 3 dimensions: [pixel, channel, waves]");
    p_dim = mxGetDimensions(M_APO);
	int P_check = (int)p_dim[0];	// number of pixels
	N_check = (int)p_dim[1];	// number of channels 
    W_check = (int)1; if (ndim == 3) W_check = (int)p_dim[2];	// number of waves
	if (P!=P_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of pixels in delay data & apodization matrix do not match.");
    if (N!=N_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of channels in channel data & apodization matrix do not match.");
    if (W!=W_check) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Number of waves in channel data & requested delays do not match.");
    if (mxIsDouble(M_APO)) mexErrMsgTxt("The apodization matrix should be single precision");    
    
    // Sampling frequency
	const size_t len_Fs = mxGetNumberOfElements(M_FS);
	if (len_Fs != 1) mexErrMsgTxt("The sampling frequency should be an escalar");
	if (mxIsDouble(M_FS)) mexErrMsgTxt("The sampling frequency should be single precision");
    const float Fs = *((float*)mxGetData(M_FS));
	const float dt = 1 / Fs;
	if (verbose) mexPrintf("Sampling frequency				%0.2f MHz\n", Fs / 1e6);

	// Initial time
	const size_t len_t0 = mxGetNumberOfElements(M_T0);
	if (len_t0 != 1) mexErrMsgTxt("The initial time should be an escalar");
	if (mxIsDouble(M_T0)) mexErrMsgTxt("The initial time should be single precision");
	const float t0 = *((float*)mxGetData(M_T0));
	if (verbose) {
        mexPrintf("Initial time					%0.2f us\n", t0*1e6);
        mexPrintf("---------------------------------------------------------------\n");
	}
    
    // modulation frequency
	float wd = 0;
	bool IQ_version = false;
	const size_t len_fd = mxGetNumberOfElements(M_FD);
	if (len_fd != 1) mexErrMsgTxt("Modulation frequency 'fd' should be an escalar");
	if (mxIsDouble(M_FD)) mexErrMsgTxt("The modulation frequency should be single precision");
    const float fd = *((float*)mxGetData(M_FD));
	if (fd > EPS) {
		if (verbose) mexPrintf("Modulation frequency:			%0.2f MHz\n", fd / 1e6);
		wd = 2 * PI*fd;
		IQ_version = true;
	}
    
    ///////////////////////////////////
	// OUTPUTS VAR
	mwSize out_size2[3];
	out_size2[0] = P;  // pixels
	out_size2[1] = F;  // frames
    out_size2[2] = W;  // waves
	M_D = mxCreateNumericArray(3, (const mwSize*)&out_size2, mxSINGLE_CLASS, mxCOMPLEX);
	float* Dr = (float*)mxGetData(M_D);
	float* Di = (float*)mxGetImagData(M_D);

    // size variables
    const unsigned int LN=L*N;
    const unsigned int LNW=LN*W;
    const unsigned int PN=P*N;
    const unsigned int PNW=PN*W;
    const unsigned int PF=P*F;
    const unsigned int PFN=PF*N;
    
	// build the 4D matrix as a vector of vector of vectors of pointers
	float* Pr = (float*)mxGetData(M_P);         // real part of channel data
	float* Pi = (float*)mxGetImagData(M_P);		// imaginary part of channel data
    std::vector<std::vector<vec_p_float>> pr;
    std::vector<std::vector<vec_p_float>> pi;
    for (int f = 0; f < F; f++) {
		std::vector<vec_p_float> temp_temp_r;
        std::vector<vec_p_float> temp_temp_i;
		for (int w = 0; w < W; w++) {
			vec_p_float temp_r; // create an array, don't work directly on buff yet.
            vec_p_float temp_i; // create an array, don't work directly on buff yet.
			for (int n = 0; n < N; n++) {
                temp_r.push_back(Pr + L*n + LN*w + LNW*f);
                temp_i.push_back(Pr + L*n + LN*w + LNW*f);
            }
			temp_temp_r.push_back(temp_r); // Store the array in the buffer
            temp_temp_i.push_back(temp_i); // Store the array in the buffer
		}
		pr.push_back(temp_temp_r); // Store the array in the buffer
        pi.push_back(temp_temp_i); // Store the array in the buffer
	}
	
    // build the delay & apodization 3D matrix as a vector of vector of pointers
    float* p_delay  = (float*)mxGetData(M_DELAY);
    float* p_apo  = (float*)mxGetData(M_APO);	
	std::vector<vec_p_float> delay;
    std::vector<vec_p_float> apo;
	for (int w = 0; w < W; w++) {
        vec_p_float temp_delay; // create an array, don't work directly on buff yet.
        vec_p_float temp_apo;   // create an array, don't work directly on buff yet.
        for (int n = 0; n < N; n++) {
            temp_delay.push_back(p_delay + P*n + PN*w);
            temp_apo.push_back(p_apo + P*n + PN*w);
        }
		delay.push_back(temp_delay); // Store the array in the buffer
        apo.push_back(temp_apo); // Store the array in the buffer
    }
        
	//////////////////////////////////////////////////////
	// Beamforming loop
	if (verbose) { mexPrintf("Beamforming started\n"); mexEvalString("drawnow;"); }
    
// #if defined (_WIN_)    
// 	Concurrency::parallel_for(0, P, [&](int pp) { // pixel loop -> WIN
// #elif defined(_UNIX_)
//     tbb::strict_ppl::parallel_for(0, P, [&](int pp) { // pixel loop -> UNIX
// #endif          
// 
//         for (int rx = 0; rx<N; rx++) for (int w = 0; w<W; w++) { // channel & wave loop

#if defined (_WIN_)    
	Concurrency::parallel_for(0, N, [&](int rx) { // channel loop -> WIN
#elif defined(_UNIX_)
    tbb::strict_ppl::parallel_for(0, N, [&](int rx) { // channel loop -> UNIX
#endif          

        for (int pp = 0; pp<P; pp++) for (int w = 0; w<W; w++) { // pixel & wave loop
    
            // references
            float& c_delay = delay[w][rx][pp];      // delay
            float& c_apo = apo[w][rx][pp];          // apodization
            
            if(c_apo>1e-3) { // we skip any apodization value smaller than -60 dB
            
                float denay = (c_delay - t0)*Fs;		// untruncated sample number 
                int n0 = (int)floor(denay);             // truncated sample number
                float b = denay - n0;                   // linear interpolation coefficient 2
                float a = 1 - b;                        // linear interpolation coefficient 1

                if (n0>0 && n0<(L - 1)) {
                    if (IQ_version) {
                        float phase = wd*c_delay;
                        float coswt = cos(phase);
                        float sinwt = sin(phase);
                        for (int f = 0; f<F; f++) { // frame loop
                            // fractional part of the delay -> linear interpolation
                            float re = a*pr[f][w][rx][n0] + b*pr[f][w][rx][n0 + 1];
                            float im = a*pi[f][w][rx][n0] + b*pi[f][w][rx][n0 + 1];
                            // apply phase change and delay
                            //Dr[pp + P*rx + PN*w + PNW*f] = c_apo*(re*coswt - im*sinwt);
                            //Di[pp + P*rx + PN*w + PNW*f] = c_apo*(im*coswt + re*sinwt);
                            Dr[pp + P*f + PF*w] += c_apo*(re*coswt - im*sinwt);
                            Di[pp + P*f + PF*w] += c_apo*(im*coswt + re*sinwt);
                            
                        }
                    } else {
                        for (int f = 0; f<F; f++) { // frame loop
                            // fractional part of the delay -> linear interpolation
                            float re = a*pr[f][w][rx][n0] + b*pr[f][w][rx][n0 + 1];
                            float im = a*pi[f][w][rx][n0] + b*pi[f][w][rx][n0 + 1];
                            // apply phase change and delay
                            //Dr[pp + P*rx + PN*w + PNW*f] = c_apo*re;
                            //Di[pp + P*rx + PN*w + PNW*f] = c_apo*im;
                            Dr[pp + P*f + PF*w] += c_apo*re;
                            Di[pp + P*f + PF*w] += c_apo*im;
                        }
                    }
                }
            }
        } // end channel & wave loop
    }); // end pixel loop 
     
	if (verbose) {
		mexPrintf("Beamforming done\n");
		mexPrintf("---------------------------------------------------\n");
		mexEvalString("drawnow;");
	}

	// freeing RAM 
	pr.clear();
	pi.clear();
    delay.clear();
    apo.clear();

	if (verbose) {
		mexPrintf("---------------------------------------------------\n");
		mexEvalString("drawnow;");
	}

	return;
	
}






