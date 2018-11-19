clear all; close all;

%filename = 'FieldII_STAI_simulated_dynamic_range.uff';
%filename = 'experimental_dynamic_range_phantom.uff';
filename = 'FieldII_STAI_dynamic_range_similar_to_exp_v4.uff';

channel_data = uff.channel_data();
channel_data.read([data_path,filesep,filename],'/channel_data');
%% Scan
scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,256).', 'z_axis', linspace(7e-3,55e-3,256).');


%%

channel_data.name = {['New ',channel_data.name{1}]}
%% Beamformer

mid = midprocess.das();
mid.channel_data=channel_data;
mid.scan=scan;
mid.dimension = dimension.transmit();
mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=1.75;
mid.transmit_apodization.window=uff.window.boxcar;
mid.transmit_apodization.f_number=1.75;

% Delay and sum on receive, then coherent compounding
b_data_tx = mid.go();

%% Calculate weights to get uniform FOV. See example
[weights,array_gain_compensation,geo_spreading_compensation] = tools.uniform_fov_weighting(mid);

%% DELAY AND SUM
mid = midprocess.das();
mid.channel_data=channel_data;
mid.scan=scan;
mid.dimension = dimension.both();
mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=1.75;
mid.transmit_apodization.window=uff.window.boxcar;
mid.transmit_apodization.f_number=1.75;

% Delay and sum on receive, then coherent compounding
b_data_das = mid.go();
%%
b_data_das.data = b_data_das.data.*weights(:); %Compensate for geometrical spreading
b_data_das.plot([],[],80);

%% COHERENCE FACTOR
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf.receive_apodization = mid.receive_apodization;
cf.transmit_apodization = mid.transmit_apodization;
cf.input = b_data_tx;
b_data_cf = cf.go();
b_data_cf.data = b_data_cf.data.*weights(:);
b_data_cf.plot();

%% MINIMUM VARIANCE
mv = postprocess.capon_minimum_variance();
mv.dimension = dimension.receive;
mv.transmit_apodization = mid.transmit_apodization;
mv.receive_apodization = mid.receive_apodization;
mv.input = b_data_tx;
mv.scan = scan;
mv.channel_data = channel_data;
mv.K_in_lambda = 1.5;
mv.L_elements = channel_data.probe.N/2;
mv.regCoef = 1/100;
b_data_mv = mv.go();
b_data_mv.data = b_data_mv.data.*weights(:);
%%
b_data_mv.plot([],[],80);

%% EBMV
ebmv=postprocess.eigenspace_based_minimum_variance();
ebmv.dimension = dimension.receive;
ebmv.input = b_data_tx;
ebmv.channel_data = channel_data;
ebmv.scan = scan;
ebmv.K_in_lambda = 1.5;
ebmv.gamma = 0.5;
ebmv.L_elements = floor(channel_data.N_elements/2);
ebmv.transmit_apodization = mid.transmit_apodization;
ebmv.receive_apodization = mid.receive_apodization;
ebmv.regCoef = 1/100;
b_data_ebmv = ebmv.go();
b_data_ebmv.data = b_data_ebmv.data.*weights(:);
b_data_ebmv.plot();
%%

diff = tools.dynamic_range_test(channel_data,b_data_das,'DAS')

%%
diff = tools.dynamic_range_test(channel_data,b_data_cf,'CF')
diff = tools.dynamic_range_test(channel_data,b_data_mv,'MV')
%%
diff =  tools.dynamic_range_test(channel_data,b_data_ebmv,'EBMV')