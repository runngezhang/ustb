clear all; close all;

filename = 'FieldII_STAI_dynamic_range_alt_2.uff';

channel_data = uff.channel_data();
channel_data.read([data_path,filesep,filename],'/channel_data');
%% Scan
%
% The scan area is defines as a collection of pixels spanning our region of
% interest. For our example here, we use the *linear_scan* structure,
% which is defined with two components: the lateral range and the
% depth range. *scan* too has a useful *plot* method it can call.

scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,1024).', 'z_axis', linspace(8e-3,55e-3,2048).');

%scan=uff.linear_scan('x_axis',linspace(-19e-3,19e-3,256).', 'z_axis', linspace(10e-3,55e-3,256).');
%scan=uff.linear_scan('x_axis',linspace(-12.5e-3,-2.5e-3,50).', 'z_axis', linspace(43e-3,46e-3,200).');
%%
% demod = preprocess.demodulation;
% demod.input = channel_data;
% 
% channel_data_demod = demod.go();
%% Beamformer

mid = midprocess.das();
mid.channel_data=channel_data;
mid.scan=scan;
mid.dimension = dimension.transmit();

mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=1.75;

mid.transmit_apodization.window=uff.window.boxcar;
mid.transmit_apodization.f_number=1.75;
%%
% Delay and sum on receive, then coherent compounding
b_data_tx = mid.go();

%% Calculate weights to get uniform FOV. See example
[weights,array_gain_compensation,geo_spreading_compensation] = tools.uniform_fov_weighting(mid);

%% Put the waits in a b_data struct to be able to save them later
b_data_weights = uff.beamformed_data();                                       
b_data_weights.scan = scan;
b_data_weights.data = weights(:);

%% Add some noise to avoid "coherent" zero-signals
b_data_tx.data = b_data_tx.data + randn(size(b_data_tx.data))*eps;

%% DELAY AND SUM
das=postprocess.coherent_compounding();
das.input = b_data_tx;
b_data_das = das.go();
das_img = b_data_das.get_image('none').*weights;  % Compensation weighting
das_img = db(abs(das_img./max(das_img(:))));                 % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
colormap gray;caxis([-120 0]);axis image;title('DAS');xlabel('x [mm]');ylabel('z [mm]');

%%
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf.receive_apodization = mid.receive_apodization;
cf.transmit_apodization = mid.transmit_apodization;
cf.input = b_data_tx;
b_data_cf = cf.go();
cf_img = b_data_cf.get_image('none').*weights;
cf_img = db(abs(cf_img./max(cf_img(:))));                 % Normalize on max
f2 = figure(2);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,cf_img);
colormap gray;caxis([-60 0]);axis image;title('CF');xlabel('x [mm]');ylabel('z [mm]');

%%
pcf = postprocess.phase_coherence_factor();
pcf.dimension = dimension.receive;
pcf.receive_apodization = mid.receive_apodization;
pcf.transmit_apodization = mid.transmit_apodization;
pcf.input = b_data_tx;
b_data_pcf = pcf.go();
pcf_img = b_data_pcf.get_image('none').*weights;
pcf_img = db(abs(pcf_img./max(pcf_img(:))));
f3 = figure(3);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,pcf_img);
colormap gray;caxis([-60 0]);axis image;title('PCF');xlabel('x [mm]');ylabel('z [mm]');

%%
gcf=postprocess.generalized_coherence_factor_OMHR();
gcf.dimension = dimension.receive;
gcf.transmit_apodization = mid.transmit_apodization;
gcf.receive_apodization = mid.receive_apodization;
gcf.input = b_data_tx;
gcf.channel_data = channel_data;
gcf.M0 = 2;
b_data_gcf = gcf.go();

gcf_img = b_data_gcf.get_image('none').*weights;
gcf_img = db(abs(gcf_img./max(gcf_img(:))));
f4 = figure(4);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img);
colormap gray;caxis([-60 0]);axis image;title('GCF');xlabel('x [mm]');ylabel('z [mm]');

%%
gcf_alt=postprocess.generalized_coherence_factor();
gcf_alt.dimension = dimension.receive;
gcf_alt.transmit_apodization = mid.transmit_apodization;
gcf_alt.receive_apodization = mid.receive_apodization;
gcf_alt.input = b_data_tx;
gcf_alt.M0 = 2;
b_data_gcf_alt = gcf_alt.go();

gcf_img_alt = b_data_gcf_alt.get_image('none').*weights;
gcf_img_alt = db(abs(gcf_img_alt./max(gcf_img_alt(:))));
f5 = figure(6);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img_alt);
colormap gray;caxis([-60 0]);axis image;title('GCF Alternative');xlabel('x [mm]');ylabel('z [mm]');

%%
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

mv_img = b_data_mv.get_image('none').*weights;
mv_img = db(abs(mv_img./max(mv_img(:))));
f5 = figure(7);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
colormap gray;caxis([-60 0]);axis image;colorbar;title('MV');xlabel('x [mm]');ylabel('z [mm]');

%% EIGENSPACE BASED MINIMUM VARIANCE
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
ebmv_img = b_data_ebmv.get_image('none').*weights;
ebmv_img = db(abs(ebmv_img./max(ebmv_img(:))));
f6 = figure(8);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,ebmv_img);
colormap gray;caxis([-60 0]);axis image;colorbar;title('EBMV');xlabel('x [mm]');ylabel('z [mm]');

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
dmas_img = b_data_dmas.get_image('none').*weights;
dmas_img = db(abs(dmas_img./max(dmas_img(:))));
f7 = figure(9);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,dmas_img);
colormap gray;caxis([-60 0]);axis image;colorbar;title('DMAS');xlabel('x [mm]');ylabel('z [mm]');

%% Gray Level Transform
glt = postprocess.gray_level_transform();
glt.a = 0.0001;
glt.b = 0;
glt.c = 1;
glt.plot_functions = 1;

glt.input = b_data_das;

b_data_glt = glt.go();
b_data_glt.plot()


%%
b_data_tx.write([data_path,filesep,filename],'/b_data_tx');
b_data_das.write([data_path,filesep,filename],'/b_data_das');
b_data_cf.write([data_path,filesep,filename],'/b_data_cf');
b_data_pcf.write([data_path,filesep,filename],'/b_data_pcf');
b_data_gcf.write([data_path,filesep,filename],'/b_data_gcf');
b_data_mv.write([data_path,filesep,filename],'/b_data_mv');
b_data_ebmv.write([data_path,filesep,filename],'/b_data_ebmv');
b_data_dmas.write([data_path,filesep,filename],'/b_data_dmas');
b_data_glt.write([data_path,filesep,filename],'/b_data_glt');
b_data_weights.write([data_path,filesep,filename],'/b_data_weights');

%%
%Using "colorbrewer" colors
colors =     [27   158   119; ...
             217    95     2; ...
             117   112   179; ...
             231    41   138; ...
             102   166    30; ...
             230   171     2; ...
             166   118    29; ...
             102   102   102]/255;
         
image.all{1} = das_img;
image.tags{1} = 'DAS';
image.all{2} = dmas_img;
image.tags{2} = 'DMAS';

[handle_2] = plotLateralLine(b_data_tx,image,45,-9,-6,colors)


%%
% %%
% addpath functions/
% 
% image.all{1} = das_img;
% image.tags{1} = 'DAS';
% image.all{5} = mv_img;
% image.tags{5} = 'MV';
% image.all{6} = ebmv_img;
% image.tags{6} = 'EBMV';
% image.all{2} = cf_img;
% image.tags{2} = 'CF';
% image.all{7} = pcf_img;
% image.tags{7} = 'PCF';
% image.all{3} = gcf_img;
% image.tags{3} = 'GCF';
% image.all{4} = das_img;%dmas_img;
% image.tags{4} = 'DMAS';
% 
% %%
% [meanLines,x_axis] = getMeanLateralLines(b_data_das,image,47.5,52.5,-20,20);
% 
% mask=abs(scan.x_axis)<10e-3;
% theory=-60*(scan.x_axis(mask)+10e-3)/30e-3-15;
% 
% f88 = figure(90);clf;hold all;
% subplot(211);
% plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1});hold on;
% plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2});
% plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3});
% plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4});
% plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
% ylim([-100 0]);
% xlim([-20 20])
% ylim([-100 0]);
% xlim([-20 20])
% legend show; grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% title('Lateral gradient');
% 
% subplot(212); hold all;
% plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5});
% plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6});
% plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7});
% plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
% ylim([-100 0]);
% xlim([-20 20])
% legend show; grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% %saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'],'eps2c')
% 
% %%
% [meanLines,x_axis] = getMeanLateralLines(b_data_das,image,47.5,52.5,-12.5,12.5);
% 
% mask=abs(scan.x_axis)<12.5e-3;
% theory=-60*(scan.x_axis(mask)+12.5e-3)/30e-3;
% 
% f88 = figure(88);clf;hold all;
% subplot(211);
% plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1});hold on;
% plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2});
% plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3});
% plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4});
% plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
% 
% ylim([-100 0]);
% xlim([-12.5 12.5])
% legend('Location','sw'); grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% title('Lateral gradient');
% 
% subplot(212);hold all;
% plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5});
% plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6});
% plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7});
% plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
% ylim([-100 0]);
% xlim([-12.5 12.5])
% legend('Location','sw'); grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% %saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'],'eps2c')
% 
% %%
% 
% [meanLines_axial,z_axis] = getMeanAxialLines(b_data_das,image,scan.z(1)*10^3,scan.z(end)*10^3,11,14);
% 
% mask= scan.z_axis<50e-3 & scan.z_axis>10e-3;
% theory=-60*(scan.z_axis(mask)-10e-3)/30e-3;
% 
% f33 = figure(33);clf; hold all;
% %set(f33,'Position',[100 100 1200 800]); hold all;
% plot(z_axis,meanLines_axial.all{1}-max(meanLines_axial.all{1}),'Linewidth',2,'DisplayName',image.tags{1});
% plot(z_axis,meanLines_axial.all{2}-max(meanLines_axial.all{2}),'Linewidth',2,'DisplayName',image.tags{2});
% plot(z_axis,meanLines_axial.all{3}-max(meanLines_axial.all{3}),'Linewidth',2,'DisplayName',image.tags{3});
% plot(z_axis,meanLines_axial.all{4}-max(meanLines_axial.all{4}),'Linewidth',2,'DisplayName',image.tags{4});
% plot(z_axis,meanLines_axial.all{5}-max(meanLines_axial.all{5}),'Linewidth',2,'DisplayName',image.tags{5});
% plot(scan.z_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
% ylim([-100 0]);
% xlim([10 43.5]);
% legend show; grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('z [mm]');
% title('Axial gradient');
% saveas(f33,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/axial_gradient'],'eps2c')
% addpath ../Functions/
% 
% %%
% 
% [meanLines,x_axis] = getMeanLateralLines(b_data_das,image,35,40,scan.x(1)*10^3,scan.x(end)*10^3);
% theoretical = [-100 ones(1,255)*0 ones(1,255)*-10 -100 ones(1,256)*-100 ones(1,256)*-100 -100 ones(1,255)*0 ones(1,255)*-25 -100];
% x_axis_theoretical = linspace(-15,0,1536);
% 
% f89 = figure(89);clf;hold all;
% subplot(211)
% plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1});hold on;
% plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2});
% plot(x_axis_theoretical,theoretical,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
% ylim([-80 0]);
% xlim([-17.5 1.5]);
% legend show; grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% title('Boxes');
% subplot(212); hold all;
% plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3});
% plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4});
% plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5});
% ylim([-80 0]);
% xlim([-17.5 1.5]);
% legend show; grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% title('Boxes');
% %saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/boxes'],'eps2c')