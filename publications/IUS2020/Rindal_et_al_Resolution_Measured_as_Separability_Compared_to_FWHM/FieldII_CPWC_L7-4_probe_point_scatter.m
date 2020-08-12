clear all; close all;

%% Channel Data
% 
% In this part of the code, we creat a uff data structure to specifically
% store the captured ultrasound channel data.

channel_data = uff.channel_data();
%%
filename = [ustb_path(),'/data/FieldII_CPWC_point_scatterers_res_v2.uff'];

%%
channel_data.read(filename,'/channel_data');

%% Scan
%
% The scan area is defines as a collection of pixels spanning our region of 
% interest. For our example here, we use the *linear_scan* structure, 
% which is defined with two components: the lateral range and the 
% depth range. *scan* too has a useful *plot* method it can call.

sca=uff.linear_scan('x_axis',linspace(-3e-3,3e-3,512).', 'z_axis', linspace(35e-3,42e-3,256).');

%% Pipeline
%
% With *channel_data* and a *scan* we have all we need to produce an
% ultrasound image. We now use a USTB structure *pipeline*, that takes an
% *apodization* structure in addition to the *channel_data* and *scan*.

pipe=pipeline();
pipe.channel_data=channel_data;
pipe.scan=sca;

pipe.receive_apodization.window=uff.window.none;
%pipe.receive_apodization.f_number=F_number;

%% 
%
% The *pipeline* structure allows you to implement different beamformers 
% by combination of multiple built-in *processes*. By changing the *process*
% chain other beamforming sequences can be implemented. It returns yet 
% another *UFF* structure: *beamformed_data*.
% 
% To achieve the goal of this example, we use delay-and-sum (implemented in 
% the *das_mex()* process) as well as coherent compounding.

b_data=pipe.go({midprocess.das()});

% Display images
b_data.plot();

%%
img = b_data.get_image();

for i = 1:7
figure;
subplot(211)
imagesc(img(:,:,i))
axis image;colormap gray;caxis([-60 0])
subplot(212)
plot(sca.x_axis*1000,img(115,:,i))
xlim([-3 3])
end
%% Save UFF dataset
% 
% Finally, we can save the data into a UFF file.

mid = midprocess.das();
mid.channel_data=channel_data;
mid.scan=sca;
mid.dimension = dimension.transmit();

mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=1.75;

mid.transmit_apodization.window=uff.window.boxcar;
mid.transmit_apodization.f_number=1.75;
%%
% Delay and sum on receive, then coherent compounding
b_data_tx = mid.go();


%% Add some noise to avoid "coherent" zero-signals
b_data_tx.data = b_data_tx.data + randn(size(b_data_tx.data))*eps;

%% DELAY AND SUM
das=postprocess.coherent_compounding();
das.input = b_data_tx;
b_data_das = das.go();
b_data_das.plot
%das_img = b_data_das.get_image('none');  % Compensation weighting
%das_img = db(abs(das_img./max(das_img(:))));                 % Normalize on max
%f1 = figure(1);clf;
%imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
%colormap gray;caxis([-120 0]);axis image;title('DAS');xlabel('x [mm]');ylabel('z [mm]');

%%
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf.receive_apodization = mid.receive_apodization;
cf.transmit_apodization = mid.transmit_apodization;
cf.input = b_data_tx;
b_data_cf = cf.go();
b_data_cf.plot
%cf_img = b_data_cf.get_image('none').*weights;
%cf_img = db(abs(cf_img./max(cf_img(:))));                 % Normalize on max
%f2 = figure(2);
%imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,cf_img);
%colormap gray;caxis([-60 0]);axis image;title('CF');xlabel('x [mm]');ylabel('z [mm]');

%%
pcf = postprocess.phase_coherence_factor();
pcf.dimension = dimension.receive;
pcf.receive_apodization = mid.receive_apodization;
pcf.transmit_apodization = mid.transmit_apodization;
pcf.input = b_data_tx;
b_data_pcf = pcf.go();
b_data_pcf.plot
%pcf_img = b_data_pcf.get_image('none').*weights;
%pcf_img = db(abs(pcf_img./max(pcf_img(:))));
%f3 = figure(3);clf;
%imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,pcf_img);
%colormap gray;caxis([-60 0]);axis image;title('PCF');xlabel('x [mm]');ylabel('z [mm]');

%%
gcf=postprocess.generalized_coherence_factor();
%gcf.dimension = dimension.receive;
%gcf.transmit_apodization = mid.transmit_apodization;
%gcf.receive_apodization = mid.receive_apodization;
gcf.input = b_data_tx;
%gcf.channel_data = channel_data;
gcf.M0 = 2;
b_data_gcf = gcf.go();
b_data_gcf.plot()
%gcf_img = b_data_gcf.get_image('none').*weights;
%gcf_img = db(abs(gcf_img./max(gcf_img(:))));
%f4 = figure(4);
%imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img);
%colormap gray;caxis([-60 0]);axis image;title('GCF');xlabel('x [mm]');ylabel('z [mm]');

%%
mv = postprocess.capon_minimum_variance();
mv.dimension = dimension.receive;
mv.transmit_apodization = mid.transmit_apodization;
mv.receive_apodization = mid.receive_apodization;
mv.input = b_data_tx;
mv.scan = sca;
mv.channel_data = channel_data;
mv.K_in_lambda = 1.5;
mv.L_elements = channel_data.probe.N/2;
mv.regCoef = 1/100;
b_data_mv = mv.go();
b_data_mv.plot()
%mv_img = b_data_mv.get_image('none').*weights;
%mv_img = db(abs(mv_img./max(mv_img(:))));
%f5 = figure(7);clf;
%imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
%colormap gray;caxis([-60 0]);axis image;colorbar;title('MV');xlabel('x [mm]');ylabel('z [mm]');

%% EIGENSPACE BASED MINIMUM VARIANCE
ebmv=postprocess.eigenspace_based_minimum_variance();
ebmv.dimension = dimension.receive;
ebmv.input = b_data_tx;
ebmv.channel_data = channel_data;
ebmv.scan = sca;
ebmv.K_in_lambda = 1.5;
ebmv.gamma = 0.5;
ebmv.L_elements = floor(channel_data.N_elements/2);
ebmv.transmit_apodization = mid.transmit_apodization;
ebmv.receive_apodization = mid.receive_apodization;
ebmv.regCoef = 1/100;

b_data_ebmv = ebmv.go();

%% Process DMAS
dmas=postprocess.delay_multiply_and_sum();
dmas.dimension = dimension.receive;
dmas.transmit_apodization = mid.transmit_apodization;
dmas.receive_apodization = mid.receive_apodization;
dmas.input = b_data_tx;
dmas.channel_data = channel_data;
b_data_dmas = dmas.go();
b_data_dmas.plot(6,['DMAS'])


%%
b_data_das.write([filename],'/b_data_das');
b_data_cf.write([filename],'/b_data_cf');
b_data_pcf.write([filename],'/b_data_pcf');
b_data_gcf.write([filename],'/b_data_gcf');
b_data_mv.write([filename],'/b_data_mv');
b_data_ebmv.write([filename],'/b_data_ebmv');
b_data_dmas.write([filename],'/b_data_dmas');