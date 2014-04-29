#define _USE_MATH_DEFINES

typedef struct
{
    float tx_angle;
    float rx_angle;
    float c0;
    float t0;
    float f_demod;
    float FN;
    float dx;
    float fs_in;
    float fs_out;    
} BF_Params;

#ifdef _USE_TRIG_TABLES_  
	const int trig_table_sz = _TRIG_TABLE_SZ_;
#endif

#ifdef _USE_V4SF_   
	void init_trig_tables_v4sf();
	void planewave_beamforming_rf_v4sf(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_data,const int ch_data_sz[2],float* rf_bf_data,BF_Params params,const float* apod_mask);

	void planewave_beamforming_iq_v4sf(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_dataI,float* ch_dataQ,const int ch_data_sz[2],float* iq_bf_dataI,float* iq_bf_dataQ,BF_Params params,const float* apod_mask);
#else
	void init_trig_tables();
	void planewave_beamforming_rf(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_data,const int ch_data_sz[2],float* rf_bf_data,BF_Params params,const float* apod_mask);

	void planewave_beamforming_iq(const float* x_coords,const int Nx,const float* r_coords,const int Nr,float* ch_dataI,float* ch_dataQ,const int ch_data_sz[2],float* iq_bf_dataI,float* iq_bf_dataQ,BF_Params params,const float* apod_mask);
#endif
