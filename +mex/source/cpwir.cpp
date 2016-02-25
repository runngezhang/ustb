#include "stdafx.h"
#include <math.h>
#include "C:\Program Files\MATLAB\R2014b\extern\include\matrix.h"
#include "C:\Program Files\MATLAB\R2014b\extern\include\mex.h"
#include <vector>
#include <concurrent_vector.h>
#include <tuple>
#include <valarray>
#include <set>
#include <ppl.h>

// compulsory input
#define	M_P			prhs[0]	// CPW dataset [samples, channels, firings, frames]
#define	M_X			prhs[1]	// vector of x components [pixels] (m)
#define	M_Z			prhs[2]	// vectors of z components [pixels] (m)
#define	M_ELE		prhs[3]	// geometry of the probe [x, z] (m, m)
#define	M_C			prhs[4] // speed of sound (m/s)
#define	M_ANGLES    prhs[5] // vector or angles [rad]
#define	M_TX_APO	prhs[6] // Transmitting apodization [pixels, angles] 
#define	M_RX_APO	prhs[7] // Receiving apodization [pixels, channels]
#define	M_FS		prhs[8] // sampling frequency (Hz)
#define M_T0		prhs[9]	// initial time (s)
#define	M_FD		prhs[10] // modulation frequency (Hz)

// optional input
#define	M_VERBOSE	prhs[11] // verbose flag

// output
#define	M_D			plhs[0] // output data [pixels, frames]

#define EPS 1e-6
#define PI 3.14159265359

// concurrent types
#define T_c_im_tuple std::tuple<int,std::vector<std::vector<std::vector<double>>>>
#define T_c_pw_tuple  std::tuple<int,std::vector<std::vector<double>>>

// types
typedef std::vector<double*> vec_p_double;

//-----------------------------------------------------------------------------

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
		
	///////////////////////////////////////////////////
	// measurement of elapsed time
	LARGE_INTEGER ticks_second;    // ticks per second
	LARGE_INTEGER start_time, stop_time;    // ticks
	QueryPerformanceFrequency(&ticks_second);

	///////////////////////////////
	// CHECKING ARGUMENTS
	///////////////////////////////////////
	// number of inputs
	if (nrhs<11) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nrhs", "Too few input arguments");
	//if (nrhs>11) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nrhs", "Too many input arguments");
	if (nlhs>1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nlhs", "Too many output arguments");

	///////////////////////////////////////
	// VERBOSE
	///////////////////////////////////////
	bool verbose = true;
	if (nrhs==12)  verbose = (*((double*)mxGetData(M_VERBOSE)))>EPS;

	if (verbose) {
		mexPrintf("---------------------------------------------------------------\n");
		mexPrintf(" Coherent Plane Wave Compunding Beamformer \n");
		mexPrintf("---------------------------------------------------------------\n");
		mexPrintf(" Vers:  0.8\n");
		mexPrintf(" Auth:  Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)\n");
		mexPrintf(" Dept:  ISB (NTNU)\n");
		mexPrintf(" Date:  19/11/2015\n");
		mexPrintf("---------------------------------------------------------------\n");
	}
	
	///////////////////////////////////////
	// IQ or RF -> modulation frequency
	///////////////////////////////////////
	double wd = 0;
	bool IQ_version=false;
	const size_t len_fd = mxGetNumberOfElements(M_FD);
	if (len_fd != 1) mexErrMsgTxt("Modulation frequency 'fd' should be an escalar");
	const double fd = *((double*)mxGetData(M_FD));
	if (fd > EPS) {
		if(verbose) mexPrintf("Modulation frequency:			%0.2f MHz\n", fd / 1e6);
		wd = 2 * PI*fd;
		IQ_version = true;
	}

	/////////////////////////////////////
	// DATASET
	///////////////////////////////////////
	int ndim = (int)mxGetNumberOfDimensions(M_P);
	if (ndim<2 || ndim>4) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown STA dataset format. Expected [samples, channels, firings, frames]");
	const mwSize* p_dim = mxGetDimensions(M_P);
	int L = (int)p_dim[0];	// number of time samples
	int N = (int)p_dim[1];	// number of channels 
	int NA = (int)1;		// number of angles 
	int F = (int)1;			// number of frames
	if (ndim == 4) {
		NA = (int)p_dim[2]; // number of angles
		F = (int)p_dim[3];	// number of frames
	}
	if (ndim == 3) {
		NA = (int)p_dim[2]; // number of angles
	}

	if (verbose) {
		mexPrintf("Samples							%i\n", L);
		mexPrintf("Channels						%i\n", N);
		mexPrintf("Plane waves						%i\n", NA);
		mexPrintf("Frames							%i\n", F);
	}

	// build the real 4D matrix as a vector of vectors of pointers
	double* Pr = mxGetPr(M_P);	// real part of STA signals
	std::vector<std::vector<vec_p_double>> pr;
	for (int f = 0; f < F; f++) {
		std::vector<vec_p_double> temp_temp_r;
		for (int i = 0; i < NA; i++) {
			vec_p_double temp_r; // create an array, don't work directly on buff yet.
			for (int j = 0; j < N; j++) temp_r.push_back(Pr + L*j + N*L*i + NA*N*L*f);
			temp_temp_r.push_back(temp_r); // Store the array in the buffer
		}
		pr.push_back(temp_temp_r); // Store the array in the buffer
	}

	// build imaginary 4D matrix as a vector of vectors of pointers
	double* Pi = mxGetPi(M_P);					// imaginary part of STA signals
	std::vector<std::vector<vec_p_double>> pi;
	if (IQ_version) for (int f = 0; f < F; f++) {
		std::vector<vec_p_double> temp_temp_i;
		for (int i = 0; i < NA; i++) {
			vec_p_double temp_i;
			for (int j = 0; j < N; j++) temp_i.push_back(Pi + L*j + N*L*i + NA*N*L*f);
			temp_temp_i.push_back(temp_i);
		}
		pi.push_back(temp_temp_i);
	}
	

	///////////////////////////////////////
	// x vector
	///////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_X);
	//mexPrintf("x dimensions %i \n",ndim);
	if (ndim != 2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "The x vector should have 2 dimensions");
	const mwSize* f_dim = mxGetDimensions(M_X);
	const int P = (int)f_dim[1];
	//mexPrintf("x length %i \n",Lx);
	if (P < 1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "The x vector should have length > 0");
	double* x = (double*)mxGetData(M_X);

	///////////////////////////////////////
	// z vector
	///////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_Z);
	//mexPrintf("z dimensions %i \n",ndim);
	if (ndim != 2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "The z vector should have 2 dimensions");
	f_dim = mxGetDimensions(M_Z);
	if (P != f_dim[1]) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "The dimensions of x and z vectors should agree");
	double* z = (double*)mxGetData(M_Z);
	//for(int i=0;i<Lz;i++) mexPrintf("z %i -> %0.2f\n",i+1,z[i]*1e3);
	if(verbose) mexPrintf("Pixels							%i\n", P);

	///////////////////////////////////////
	// angles vector
	///////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_ANGLES);
	if (ndim != 2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "The angles vector should have 2 dimensions");
	f_dim = mxGetDimensions(M_ANGLES);
	if (NA != f_dim[1]) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "The length of the angles vector should agree with the dataset size");
	double* ang = (double*)mxGetData(M_ANGLES);

	///////////////////////////////////////
	// elements
	///////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_ELE);
	if (ndim != 2) mexErrMsgTxt("Probe geometry must have 2 dimensions");
	const mwSize* ele_dim = mxGetDimensions(M_ELE);
	if (ele_dim[1] != 2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Wrong dimensions of probe geometry matrix. [x y] coordinates are expected.");
	if (ele_dim[0] != N) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Wrong dimensions of probe geometry matrix. It does not agree with dataset.");
	double* ele = (double*)mxGetData(M_ELE);
	double* x0 = ele;
	double* z0 = ele + N;

	///////////////////////////////////////
	// Apodization
	/////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_TX_APO);
	if (ndim != 2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown transmit apodization format. Expected [pixels, channels]");
	const mwSize* tx_apo_dim = mxGetDimensions(M_TX_APO);
	if (tx_apo_dim[0] != P) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Transmit apodization pixels should agree, %i != %i\n", tx_apo_dim[0], P);
	if (tx_apo_dim[1] != NA) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Transmit apodization angles should agree, %i != %i\n", tx_apo_dim[1], NA);
	double* P_an_apo = (double*)mxGetData(M_TX_APO);

	ndim = (int)mxGetNumberOfDimensions(M_RX_APO);
	if (ndim != 2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Unknown receive apodization format. Expected [pixels, channels]");
	const mwSize* rx_apo_dim = mxGetDimensions(M_RX_APO);
	if (rx_apo_dim[0] != P) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Receive apodization pixels should agree, %i != %i\n", rx_apo_dim[0], P);
	if (rx_apo_dim[1] != N) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions", "Receive apodization channels should agree, %i != %i\n", rx_apo_dim[1], N);
	double* P_rx_apo = (double*)mxGetData(M_RX_APO);

	// build the transmit and receive apodization 2D matrices as vectors of pointers
	vec_p_double an_apo;
	vec_p_double rx_apo;
	for (int j = 0; j < NA; j++) an_apo.push_back(P_an_apo + P*j);
	for (int j = 0; j < N; j++) rx_apo.push_back(P_rx_apo + P*j);

	///////////////////////////////////////
	// SINGLE PARAMETERS
	///////////////////////////////////////

	// Sampling frequency
	const size_t len_Fs = mxGetNumberOfElements(M_FS);
	if (len_Fs != 1) mexErrMsgTxt("The sampling frequency should be an escalar");
	const double Fs = *((double*)mxGetData(M_FS));
	const double dt = 1 / Fs;
	if(verbose) mexPrintf("Sampling frequency				%0.2f MHz\n", Fs / 1e6);

	// Initial time
	const size_t len_t0 = mxGetNumberOfElements(M_T0);
	if (len_t0 != 1) mexErrMsgTxt("The initial time should be an escalar");
	const double t0 = *((double*)mxGetData(M_T0));
	if(verbose) mexPrintf("Initial time					%0.2f us\n", t0*1e6);

	// speed of sound
	const size_t len_c = mxGetNumberOfElements(M_C);
	if (len_c != 1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:size", "Sound speed 'c' should be an escalar");
	const double c0 = *((double*)mxGetData(M_C));
	if(verbose) mexPrintf("Speed of sound:					%0.2f m/s\n", c0);

	///////////////////////////////////////////////////////
	//// Creating data structures
	if (verbose) {
		mexPrintf("1.- Creating concurrent data structures\n"); mexEvalString("drawnow;");
		QueryPerformanceCounter(&start_time);
	}

	// creating concurrent structures for real data
	Concurrency::concurrent_vector<std::vector<double>> c_im_r;
	for (int f = 0; f < F; f++) {
		std::vector<double> im_r(P, 0.0); // vector of doubles
		c_im_r.push_back(im_r);
	}
	// creating concurrent structures for imag data
	Concurrency::concurrent_vector<std::vector<double>> c_im_i;
	if (IQ_version) for (int f = 0; f < F; f++) {
		std::vector<double> im_i(P, 0.0);
		c_im_i.push_back(im_i);
	}

	//////////////////////////////////////////////////////
	// Beamforming loop
	if (verbose) { mexPrintf("2.- Beamforming (multi-CPU)\n"); mexEvalString("drawnow;"); }
	Concurrency::parallel_for(0, P, [&](int pp) {
		// for(int nz=0; nz<Lz; nz++) { // z vector
		double& zz = z[pp];
		double& xx = x[pp];

		for (int an = 0; an < NA; an++) { // na plane wave
			double Da = (zz*cos(ang[an]) + xx*sin(ang[an])) / c0;

			for (int rx = 0; rx<N; rx++) { // rx channel
				double apo = an_apo[an][pp] * rx_apo[rx][pp];
				if (apo>1e-3) {
					double xx2 = xx - x0[rx]; xx2 *= xx2;
					double zz2 = zz - z0[rx]; zz2 *= zz2;
					double delay = Da + sqrt(xx2 + zz2) / c0;

					double denay = (delay - t0)*Fs;		// untruncated sample number 
					int n0 = (int)floor(denay);		// truncated sample number
					double b = denay - n0;				// linear interpolation coefficient 2
					double a = 1 - b;					// linear interpolation coefficient 1

					if (n0 > 0 && n0 < (L - 1)) {
						if (IQ_version) {
							double phase = wd*delay;
							double coswt = cos(phase);
							double sinwt = sin(phase);
							for (int f = 0; f < F; f++) { // frame vector
								// fractional part of the delay -> linear interpolation
								double re = a*pr[f][an][rx][n0] + b*pr[f][an][rx][n0 + 1];
								double im = a*pi[f][an][rx][n0] + b*pi[f][an][rx][n0 + 1];
								// apply phase change and delay
								c_im_r[f][pp] += apo*(re*coswt - im*sinwt);
								c_im_i[f][pp] += apo*(im*coswt + re*sinwt);
							}
						}
						else {
							for (int f = 0; f < F; f++) { // frame vector
								double re = a*pr[f][an][rx][n0] + b*pr[f][an][rx][n0 + 1];
								c_im_r[f][pp] += apo*re;
							}
						}
					}
				}
			} // end rx loop
		} // end tx loop
	}); // end pixel loop <- concurrent
	if(verbose) QueryPerformanceCounter(&stop_time);

	if (verbose) {
		mexPrintf("3.- Copying results to output structures\n");
		mexPrintf("---------------------------------------------------\n");
		mexEvalString("drawnow;");
	}

	// freeing RAM as soon as possible
	pr.clear();
	pi.clear();

	///////////////////////////////////
	// OUTPUTS VAR
	mwSize out_size2[2];
	out_size2[0] = P; // samples
	out_size2[1] = F;  // frames
	M_D = mxCreateNumericArray(2, (const mwSize*)&out_size2, mxDOUBLE_CLASS, mxCOMPLEX);
	double* Dr = mxGetPr(M_D);
	double* Di = mxGetPi(M_D);

	////////////////////////////////////////
	// transfer results to output matrices
	for (int f = 0; f<F; f++) for (int pp = 0; pp<P; pp++) Dr[pp + f*P] = c_im_r[f][pp];
	c_im_r.clear();

	if (IQ_version) for (int f = 0; f<F; f++) for (int pp = 0; pp<P; pp++) Di[pp + f*P] = c_im_i[f][pp];
	c_im_i.clear();

	if (verbose) {
		mexPrintf("Done in %0.2fs\n", (float)(stop_time.QuadPart - start_time.QuadPart) / (float)ticks_second.QuadPart);
		mexPrintf("---------------------------------------------------\n");
		mexEvalString("drawnow;");
	}

	return;
}






