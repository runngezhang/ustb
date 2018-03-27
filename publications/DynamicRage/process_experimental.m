clear all;
close all;

%%
filename = [data_path,filesep,'experimental_gradient_phantom.uff'];
channel_data = uff.channel_data();
channel_data.read(filename,'/channel_data')

%%
%%
scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,512*2).','z_axis',linspace(8e-3,52e-3,2048).')
%scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,256).','z_axis',linspace(8e-3,52e-3,256).')

%% Beamformer
%
% With *channel_data* and a *scan* we have all we need to produce an
% ultrasound image. We now use a USTB structure *beamformer*, that takes an
% *apodization* structure in addition to the *channel_data* and *scan*.

mid = midprocess.das();
mid.channel_data=channel_data;
mid.scan=scan;
mid.dimension = dimension.transmit();

mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=1.75;

mid.transmit_apodization.window=uff.window.boxcar;
mid.transmit_apodization.f_number=1.75;


%%
[weights,array_gain_compensation,geo_spreading_compensation] = ...
                                           tools.uniform_fov_weighting(mid);
%% 
b_data_weights = uff.beamformed_data();                                       
b_data_weights.scan = scan;
b_data_weights.data = weights(:);
%%

% Delay and sum on receive, then coherent compounding
b_data_tx = mid.go();


%% DELAY AND SUM
das=postprocess.coherent_compounding();
das.input = b_data_tx;
b_data_das = das.go();
f1 = figure(3);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_das.get_image);
colormap gray;caxis([-60 0]);axis image;title('DAS');xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);

%%
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf.receive_apodization = mid.receive_apodization;
cf.transmit_apodization = mid.transmit_apodization;
cf.input = b_data_tx;
b_data_cf = cf.go();
f2 = figure(2);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_cf.get_image);
colormap gray;caxis([-60 0]);axis image;title('CF');xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);

%%
pcf = postprocess.phase_coherence_factor();
pcf.dimension = dimension.receive;
pcf.receive_apodization = mid.receive_apodization;
pcf.transmit_apodization = mid.transmit_apodization;
pcf.input = b_data_tx;
b_data_pcf = pcf.go();
f3 = figure(3);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_pcf.get_image);
colormap gray;caxis([-60 0]);axis image;title('PCF');xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);

%%
gcf=postprocess.generalized_coherence_factor_OMHR();
gcf.dimension = dimension.receive;
gcf.transmit_apodization = mid.transmit_apodization;
gcf.receive_apodization = mid.receive_apodization;
gcf.input = b_data_tx;
gcf.channel_data = channel_data;
gcf.M0 = 2;
b_data_gcf = gcf.go();
f4 = figure(4);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_gcf.get_image);
colormap gray;caxis([-60 0]);axis image;title('GCF');xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);

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
%%
f5 = figure(6);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_mv.get_image);
colormap gray;caxis([-60 0]);axis image;title('MV');xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);

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
dmas_img = b_data_dmas.get_image('none');
dmas_img = db(abs(dmas_img./max(dmas_img(:))));
f7 = figure(7);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,dmas_img);
colormap gray;caxis([-60 0]);axis image;title('DMAS');xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);

%% EIGENSPACE BASED MINIMUM VARIANCE
ebmv=postprocess.eigenspace_based_minimum_variance();
ebmv.dimension = dimension.receive;
ebmv.input = b_data_tx;
ebmv.channel_data = channel_data;
ebmv.scan = scan;
ebmv.K_in_lambda = 1.5;
ebmv.gamma = 0.5;
ebmv.L_elements = floor(channel_data.probe.N/2);
ebmv.transmit_apodization = mid.transmit_apodization;
ebmv.receive_apodization = mid.receive_apodization;
ebmv.regCoef = 1/100;

b_data_ebmv = ebmv.go();
%%

%b_data_ebmv.plot(f6,['EBMV'])
ebmv_img = b_data_ebmv.get_image('none');
ebmv_img = db(abs(ebmv_img./max(ebmv_img(:))));
f6 = figure(8);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,ebmv_img);
colormap gray;caxis([-60 0]);axis image;title('EBMV');xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);


%% Gray Level Transform
glt = postprocess.gray_level_transform();
glt.input = b_data_das;

b_data_glt = glt.go();
b_data_glt.plot()

%%

f100 = figure(100);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gamma_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f100,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/GAMMA'],'eps2c')

%%
addpath functions/

b_data_tx.write(filename,'/b_data_tx');
b_data_das.write(filename,'/b_data_das');
b_data_cf.write(filename,'/b_data_cf');
b_data_pcf.write(filename,'/b_data_pcf');
b_data_gcf.write(filename,'/b_data_gcf');
b_data_mv.write(filename,'/b_data_mv');
b_data_ebmv.write(filename,'/b_data_ebmv');
b_data_dmas.write(filename,'/b_data_dmas');
b_data_glt.write(filename,'/b_data_glt');
%%

b_data_weights.write(filename,'/b_data_weights');
%%
image.all{1} = b_data_das.get_image;
image.tags{1} = 'DAS';
image.all{5} = b_data_mv.get_image;
image.tags{5} = 'MV';
image.all{6} = b_data_ebmv.get_image
image.tags{6} = 'EBMV';
image.all{2} = b_data_cf.get_image
image.tags{2} = 'CF';
image.all{7} = b_data_pcf.get_image
image.tags{7} = 'PCF';
image.all{3} = b_data_gcf.get_image
image.tags{3} = 'GCF';
image.all{4} =  b_data_das.get_image;%b_data_dmas.get_image
image.tags{4} = 'DMAS';
%%
addpath functions
addpath cbrewer

load colorbrewer.mat

colors = colorbrewer.qual.Dark2{8}./255;
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,40,50,-15,15);

mask=abs(scan.x_axis)<15e-3;
theory=-50*(scan.x_axis(mask)+15e-3)/30e-3;

f89 = figure(89);clf
subplot(211);hold all;
plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',colors(5,:));
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));

ylim([-60 0]);
xlim([-15 15])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
title('Lateral gradient');
set(gca,'Fontsize',14)

subplot(212);hold all;
plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',colors(5,:));
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(6,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(7,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(8,:));
ylim([-60 0]);
xlim([-15 15])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
title('Lateral gradient');

set(gca,'Fontsize',14)
%saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/gradient'],'eps2c')

%%
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,40,50,-20,20);
mask=abs(scan.x_axis)<15e-3;
theory=-50*(scan.x_axis(mask)+15e-3)/30e-3;

f89 = figure(89);clf;hold all;
plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',colors(5,:));
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(6,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(7,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(8,:));
ylim([-80 0]);
xlim([-20 20])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
title('Lateral gradient');

set(gca,'Fontsize',20)
%saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/gradient_2'],'eps2c')

%%
mask=abs(scan.x_axis)<15e-3;
%ix=mean(im,1);
%px=polyfit(gradient_scan.x_axis(mask),ix(mask).',1);
theory=-10-50*(scan.x_axis(mask)+15e-3)/30e-3;

%figure;

%% check gradient
gradient_scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,512).','z_axis',linspace(40e-3,50e-3,128).');
mid.scan = gradient_scan;
b_data = mid.go(); 
b_data.plot([],'B-mode',50);

im=b_data.get_image();
mask=abs(gradient_scan.x_axis)<15e-3;
ix=mean(im,1);
px=polyfit(gradient_scan.x_axis(mask),ix(mask).',1);
theory=-10-50*(gradient_scan.x_axis(mask)+15e-3)/30e-3;

figure;
plot(gradient_scan.x_axis*1e3,ix); grid on; hold on;
plot(gradient_scan.x_axis(mask)*1e3,theory,'k-'); 
plot(gradient_scan.x_axis(mask)*1e3,gradient_scan.x_axis(mask)*px(1)+px(2),'r--');
xlabel('x [mm]');
ylabel('Signal Intensity [dB]');

