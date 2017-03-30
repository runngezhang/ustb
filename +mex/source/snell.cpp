#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>
#include <concurrent_vector.h>
#include <tuple>
#include <valarray>
#include <set>
#if defined(_WIN_)
    #include <ppl.h>            // Requires VS2010+
#elif defined(_UNIX_)
    #include <parallel_for.h>   // Requires Intel tbb
#endif

// compulsory input
#define	M_P			prhs[0]	// STA dataset [samples, channels, firings, frames]
#define	M_X			prhs[1]	// vector of x components [pixels] (m)
#define	M_Z			prhs[2]	// vectors of z components [pixels] (m)
#define	M_ELE		prhs[3]	// geometry of the probe [x, z] (m, m)
#define	M_C			prhs[4] // speed of sound (m/s)
#define	M_TX_APO	prhs[5] // Transmitting apodization [pixels, channels] 
#define	M_RX_APO	prhs[6] // Receiving apodization [pixels, channels]
#define	M_FS		prhs[7] // sampling frequency (Hz)
#define M_T0		prhs[8]	// initial time (s)
#define	M_FD		prhs[9]	// modulation frequency (Hz)
#define	M_FNUMBER	prhs[10]	// F-number

// output
#define	M_D			plhs[0] // output data [pixels, snell, frames]

// constants
#define EPS 1e-6
#define PI 3.14159265359
#define PI2 3.14159265359/2

// concurrent types
#define T_c_im_tuple std::tuple<int,std::vector<std::vector<std::vector<double>>>>
#define T_c_pw_tuple  std::tuple<int,std::vector<std::vector<double>>>

//types
typedef std::vector<double*> vec_p_double;

//-----------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) { 
	
	mexPrintf("---------------------------------------------------------------\n");
	mexPrintf(" snell (STAI snell beamformer) \n");
	mexPrintf("---------------------------------------------------------------\n");
	mexPrintf(" Vers:  0.11\n");
	mexPrintf(" Auth:  Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)\n");
	mexPrintf(" Dept:  ISB (NTNU)\n");
	mexPrintf(" Date:  31/10/2016\n");
	mexPrintf("---------------------------------------------------------------\n");

	///////////////////////////////
    // CHECKING ARGUMENTS
	///////////////////////////////////////
	// number of inputs
	if (nrhs<9) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nrhs","Too few input arguments"); 
	if (nrhs>11) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nrhs","Too many input arguments"); 
	if (nlhs>1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:nlhs", "Too many output arguments"); 

	bool IQ_version=false;
	double F_number=1;

	if(nrhs==9) {
		IQ_version = false;
		mexPrintf("Data Format						RF\n");
	} else {
        const size_t len_fd = mxGetNumberOfElements(M_FD);
		if (len_fd!=1) mexErrMsgTxt("Modulation frequency 'fd' should be an escalar");
		const double fd=*((double*)mxGetData(M_FD)); 
        if (fd>EPS) {
    		IQ_version = true;
            mexPrintf("Data Format						IQ\n");
        } else {
            IQ_version = false;
            mexPrintf("Data Format						RF\n");
        }
	}

	if(nrhs==11) {
	}

	/////////////////////////////////////
	// STA MATRIX 
	///////////////////////////////////////
	int ndim = (int)mxGetNumberOfDimensions(M_P);
	if(ndim<3||ndim>4) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Unknown STA dataset format. Expected [samples, channels, firings, frames]"); 
	const mwSize* p_dim=mxGetDimensions(M_P);
	int L=(int)p_dim[0]; // number of time samples
	int N=(int)p_dim[1]; // number of channels 
	if(p_dim[2]!=N) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","In STA datasets channels and firing must agree, %i != %i\n",N,p_dim[2]); 
	int F=(int)1;
	if (ndim==4) F=(int)p_dim[3]; // number of frames
	
	mexPrintf("Samples							%i\n",L);
	mexPrintf("Channels						%i\n",N);
	mexPrintf("Firings							%i\n",N);
	mexPrintf("Frames							%i\n",F);
	
	// build the real 4D matrix as a vector of vectors of pointers
    double* Pr=mxGetPr(M_P);	// real part of STA signals
	std::vector<std::vector<vec_p_double>> pr;
	for(int f = 0; f < F; f++) {
		std::vector<vec_p_double> temp_temp_r;
		for(int i = 0; i < N; i++) {
			vec_p_double temp_r; // create an array, don't work directly on buff yet.
			for(int j = 0; j < N; j++) temp_r.push_back(Pr+L*j+N*L*i+N*N*L*f); 
			temp_temp_r.push_back(temp_r); // Store the array in the buffer
		}
		pr.push_back(temp_temp_r); // Store the array in the buffer
	}
	
	// build imaginary 4D matrix as a vector of vectors of pointers
	double* Pi=mxGetPi(M_P);					// imaginary part of STA signals
	std::vector<std::vector<vec_p_double>> pi;
	if (IQ_version) for(int f = 0; f < F; f++) {
		std::vector<vec_p_double> temp_temp_i;
		for(int i = 0; i < N; i++) {
			vec_p_double temp_i; 
			for(int j = 0; j < N; j++) temp_i.push_back(Pi+L*j+N*L*i+N*N*L*f); 
			temp_temp_i.push_back(temp_i); 
		}
		pi.push_back(temp_temp_i); 
	}

	///////////////////////////////////////
	// x vector
	///////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_X);
	//mexPrintf("x dimensions %i \n",ndim);
	if(ndim!=2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","The x vector should have 2 dimmension"); 
	const mwSize* f_dim=mxGetDimensions(M_X);
	const int P=(int)f_dim[1];
	//mexPrintf("x length %i \n",Lx);
	if(P<1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","The x vector should have length > 0"); 
	double* x=(double*)mxGetData(M_X); 
	
	///////////////////////////////////////
	// z vector
	///////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_Z);
	//mexPrintf("z dimensions %i \n",ndim);
	if(ndim!=2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","The z vector should have 2 dimmension"); 
	f_dim=mxGetDimensions(M_Z);
	if(P!=f_dim[1]) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","The dimmensions of x and z vectors should agree"); 
	double* z=(double*)mxGetData(M_Z); 
	//for(int i=0;i<Lz;i++) mexPrintf("z %i -> %0.2f\n",i+1,z[i]*1e3);
	mexPrintf("Pixels							%i\n",P);

	///////////////////////////////////////
	// elements
	///////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_ELE);
	if(ndim!=2) mexErrMsgTxt("Probe geometry must have 2 dimensions"); 
	const mwSize* ele_dim=mxGetDimensions(M_ELE);
	if(ele_dim[1]!=2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Wrong dimensions of probe geometry matrix. [x y] coordinates are expected."); 
	if(ele_dim[0]!=N) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Wrong dimensions of probe geometry matrix. It does not agree with dataset."); 
	double* ele=(double*)mxGetData(M_ELE); 
	double* x0=ele;
	double* z0=ele+N;

	///////////////////////////////////////
	// Apodization
	/////////////////////////////////////
	ndim = (int)mxGetNumberOfDimensions(M_TX_APO);
	if(ndim!=2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Unknown transmit apodization format. Expected [pixels, channels]"); 
	const mwSize* tx_apo_dim=mxGetDimensions(M_TX_APO);
	if(tx_apo_dim[0]!=P) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Transmit apodization pixels should agree, %i != %i\n",tx_apo_dim[0],P); 
	if(tx_apo_dim[1]!=N) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Transmit apodization channels should agree, %i != %i\n",tx_apo_dim[1],N); 
	double* P_tx_apo=(double*)mxGetData(M_TX_APO);	

	ndim = (int)mxGetNumberOfDimensions(M_RX_APO);
	if(ndim!=2) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Unknown receive apodization format. Expected [pixels, channels]"); 
	const mwSize* rx_apo_dim=mxGetDimensions(M_RX_APO);
	if(rx_apo_dim[0]!=P) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Receive apodization pixels should agree, %i != %i\n",rx_apo_dim[0],P); 
	if(rx_apo_dim[1]!=N) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:Dimensions","Receive apodization channels should agree, %i != %i\n",rx_apo_dim[1],N); 
	double* P_rx_apo=(double*)mxGetData(M_RX_APO);	
		
	// build the transmit and receive apodization 2D matrices as vectors of pointers
    vec_p_double tx_apo;
	vec_p_double rx_apo;
	for(int j = 0; j < N; j++) {
		tx_apo.push_back(P_tx_apo+P*j); 
		rx_apo.push_back(P_rx_apo+P*j); 
	}

	///////////////////////////////////////
	// SINGLE PARAMETERS
	///////////////////////////////////////
	// modulation frequecy
	double wd=0;
	if (IQ_version) {
		const size_t len_fd = mxGetNumberOfElements(M_FD);
		if (len_fd!=1) mexErrMsgTxt("Modulation frequency 'fd' should be an escalar");
		const double fd=*((double*)mxGetData(M_FD)); 
		mexPrintf(			"Modulation frequency:			%0.2f MHz\n", fd/1e6);
		wd=2*PI*fd;
        if(wd<1e-1) IQ_version=1;
	}

	// Sampling frequency
	const size_t len_Fs = mxGetNumberOfElements(M_FS);
	if (len_Fs!=1) mexErrMsgTxt("The sampling frequency should be an escalar");
	const double Fs=*((double*)mxGetData(M_FS)); 
	const double dt=1/Fs;
	mexPrintf(			"Sampling frequency				%0.2f MHz\n", Fs/1e6);
	
	// Initial time
	const size_t len_t0 = mxGetNumberOfElements(M_T0);
	if (len_t0!=1) mexErrMsgTxt("The initial time should be an escalar");
	const double t0=*((double*)mxGetData(M_T0)); 
	mexPrintf(			"Initial time					%0.2f us\n", t0*1e6);
	
	// speed of sound
	const size_t len_c = mxGetNumberOfElements(M_C);
	if (len_c!=1) mexErrMsgIdAndTxt("Toolbox:SRP_SRC:size","Sound speed 'c' should be an escalar");
	const double c0=*((double*)mxGetData(M_C)); 
	mexPrintf(			"Speed of sound:					%0.2f m/s\n", c0);
	
	/////////////////////////////////////////////////////
	// Creating data structures
	mexPrintf(			"1.- Creating concurrent data structures\n"); mexEvalString("drawnow;");
	
	// creating concurrent structures for real data
    int Gmax=200;
#if defined(_WIN_)
    Concurrency::concurrent_vector<std::vector<std::vector<double>>> c_im_r;
#elif defined(_UNIX_)
    tbb::concurrent_vector<std::vector<std::vector<double>>> c_im_r;    
#endif
	for(int f=0;f<F;f++) {
		std::vector<std::vector<double>> temp_im_r; // vector of vector of doubles
		for(int g=0;g<Gmax;g++) {
			std::vector<double> im_r(P,0.0); // vector of doubles
			temp_im_r.push_back(im_r);
		}
		c_im_r.push_back(temp_im_r);
	}

	// creating concurrent structures for imag data
#if defined(_WIN_)
	Concurrency::concurrent_vector<std::vector<std::vector<double>>> c_im_i;
#elif defined(_UNIX_)
    tbb::concurrent_vector<std::vector<std::vector<double>>> c_im_i;
#endif
	if(IQ_version) for(int f=0;f<F;f++) {
		std::vector<std::vector<double>> temp_im_i;
		for(int g=0;g<Gmax;g++) {
			std::vector<double> im_i(P,0.0);
			temp_im_i.push_back(im_i);
		}
		c_im_i.push_back(temp_im_i);
	}

	//////////////////////////////////////////////////////
	// Concurrent loop
	mexPrintf("2.- Launching Multi-CPU loop\n"); mexEvalString("drawnow;");
	//Concurrency::critical_section cs;
	
	double max_angle=1.0/2.0/F_number;

	mexPrintf("max_angle:			%0.2f\n", max_angle);
     
	////////////////////////////////////////////
	// imaging loop
#if defined(_WIN_)
	Concurrency::parallel_for (0, P, [&](int pp) { 
#elif defined(_UNIX_)
     tbb::parallel_for (0, P, [&](int pp) { 
#endif
	// for(int nz=0; nz<Lz; nz++) { // z vector
		double& zz=z[pp];
		double& xx=x[pp];

		for(int tx=0;tx<N;tx++) { // tx channel
			// computing tx delay
			double xx2=xx-x0[tx]; xx2*=xx2;
			double zz2=zz; zz2*=zz;
			double tx_delay=sqrt(xx2+zz2)/c0;

			// computing transmit angle
			double alpha=atan2(xx-x0[tx],zz); // <- this is probably slow
				
			//if (abs(alpha)<max_angle) 
			{ 
				for(int rx=0;rx<N;rx++) { // rx channel
					double apo=tx_apo[tx][pp]*rx_apo[rx][pp];
					if (apo>1e-3) {
						// computing rx delay
						double xx2=xx-x0[rx]; xx2*=xx2;
						double zz2=zz; zz2*=zz;
						double rx_delay=sqrt(xx2+zz2)/c0;

						// computing receive angle
						double beta=-atan2(xx-x0[rx],zz); // <- this is probably slow

						// total delay and interpolation constants
						double delay=tx_delay+rx_delay;
						double denay=(delay-t0)*Fs;		// untruncated sample number 
						int n0=(int)floor(denay);		// truncated sample number
						double b=denay-n0;				// linear interpolation coefficient 2
						double a=1-b;					// linear interpolation coefficient 1

						// computing gamma value and interpolation constants
						double gamma=(alpha-beta)/2.0;
						double ganna=(gamma+PI2)*(float)(Gmax-1)/PI;
						int g0=(int)floor(ganna);
						double gb=ganna-g0;
						double ga=1.0-gb;

						//if (abs(beta)<max_angle) 
						{
							if(n0>0&&n0<(L-1)&&g0>-1&&g0<(Gmax-1)) {
								if(IQ_version) {
									double phase=wd*delay;
									double coswt=cos(phase);
									double sinwt=sin(phase);
									for(int f=0;f<F;f++) { // frame vector
										// fractional part of the delay -> linear interpolation
										double re=a*pr[f][tx][rx][n0]+b*pr[f][tx][rx][n0+1];
										double im=a*pi[f][tx][rx][n0]+b*pi[f][tx][rx][n0+1];
										// apply phase change and delay
										c_im_r[f][g0][pp]+=ga*apo*(re*coswt-im*sinwt);
										c_im_i[f][g0][pp]+=ga*apo*(im*coswt+re*sinwt);
										c_im_r[f][g0+1][pp]+=gb*apo*(re*coswt-im*sinwt);
										c_im_i[f][g0+1][pp]+=gb*apo*(im*coswt+re*sinwt);
									}
								} else {
									for(int f=0;f<F;f++) { // frame vector
										double re=a*pr[f][tx][rx][n0]+b*pr[f][tx][rx][n0+1]; 
										c_im_r[f][g0][pp]+=ga*apo*re;
										c_im_r[f][g0+1][pp]+=gb*apo*re;
									}
								}
							} // check limits of gamma
						} // check on RX angle
					} // apodization check
				} // end rx loop
			} // check on TX angle
		} // end tx loop
	} ); // end pixel loop <- concurrent
	     
	mexPrintf("3.- Copying results to output structures\n");
	mexPrintf("---------------------------------------------------\n");
	mexEvalString("drawnow;");

	// freeing RAM as soon as possible
	pr.clear();
	pi.clear();

	///////////////////////////////////
	// OUTPUTS VAR
	mwSize out_size2[3];
	out_size2[0]= P; // samples
	out_size2[1]= Gmax; // gamma
	out_size2[2]= F;  // frames
	M_D=mxCreateNumericArray(3, (const mwSize*) &out_size2, mxDOUBLE_CLASS, mxCOMPLEX);
	double* Dr=mxGetPr(M_D);
	double* Di=mxGetPi(M_D);

	////////////////////////////////////////
	// transfer results to output matrices
	for(int f=0;f<F;f++) for(int g=0;g<Gmax;g++) for(int pp=0;pp<P;pp++) Dr[pp+g*P+f*Gmax*P]=c_im_r[f][g][pp];
	c_im_r.clear();

	if(IQ_version) for(int f=0;f<F;f++) for(int g=0;g<Gmax;g++) for(int pp=0;pp<P;pp++) Di[pp+g*P+f*Gmax*P]=c_im_i[f][g][pp];
	c_im_i.clear();

	mexPrintf("---------------------------------------------------\n");
	mexEvalString("drawnow;");
	
    return;
}
