clear all;
close all;


filename = 'data/experimental_gradient_phantom.uff';

channel_data = uff.channel_data();
channel_data.read(filename,'/channel_data')

% demod
demod = preprocess.demodulation();
demod.input = channel_data;
channel_data_demod = demod.go();
%%
scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,512).','z_axis',linspace(8e-3,52e-3,512).');

%% Beamformer
%
% With *channel_data* and a *scan* we have all we need to produce an
% ultrasound image. We now use a USTB structure *beamformer*, that takes an
% *apodization* structure in addition to the *channel_data* and *scan*.

mid = midprocess.das();
mid.channel_data=channel_data_demod;
mid.scan=scan;
mid.dimension = dimension.transmit();

mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=2;

mid.transmit_apodization.window=uff.window.hanning;
mid.transmit_apodization.f_number=2;

% Delay and sum on receive, then coherent compounding
b_data_tx = mid.go();

%% Add some noise to avoid white areas in e.g. CF

%b_data_tx.data = b_data_tx.data + randn(size(b_data_tx.data))*eps;

%%

figure;hold all;
%plot(data_test(:,64));
plot(real(b_data_tx.data(:,64)));
%% DELAY AND SUM
das=postprocess.coherent_compounding();
das.input = b_data_tx;
b_data_das = das.go();
f1 = figure(1);clf;
%b_data_das.plot(f1,['DAS'],60)
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_das.get_image);
colormap gray;caxis([-60 0]);axis image;colorbar;title('DAS');xlabel('x [mm]');ylabel('z [mm]');
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DAS'],'eps2c')

%%
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf.receive_apodization = mid.receive_apodization;
cf.transmit_apodization = mid.transmit_apodization;
cf.input = b_data_tx;
b_data_cf = cf.go();
f2 = figure(2);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_cf.get_image);
colormap gray;caxis([-60 0]);axis image;colorbar;title('CF');xlabel('x [mm]');ylabel('z [mm]');
saveas(f2,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/CF'],'eps2c')
%%
pcf = postprocess.phase_coherence_factor();
pcf.dimension = dimension.receive;
pcf.receive_apodization = mid.receive_apodization;
pcf.transmit_apodization = mid.transmit_apodization;
pcf.input = b_data_tx;
b_data_pcf = pcf.go();
f3 = figure(3);clf;
%b_data_pcf.plot(f3,'PCF');
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_pcf.get_image);
colormap gray;caxis([-60 0]);axis image;colorbar;title('PCF');xlabel('x [mm]');ylabel('z [mm]');
saveas(f3,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/PCF'],'eps2c')

%%
gcf=postprocess.generalized_coherence_factor();
gcf.dimension = dimension.receive;
gcf.transmit_apodization = mid.transmit_apodization;
gcf.receive_apodization = mid.receive_apodization;
gcf.input = b_data_tx;
gcf.channel_data = channel_data;
gcf.nBeams = 5;
b_data_gcf = gcf.go();
f4 = figure(4);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_gcf.get_image);
colormap gray;caxis([-60 0]);axis image;colorbar;title('GCF');xlabel('x [mm]');ylabel('z [mm]');
saveas(f4,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/GCF'],'eps2c')

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
colormap gray;caxis([-60 0]);axis image;colorbar;title('MV');xlabel('x [mm]');ylabel('z [mm]');
saveas(f5,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/MV'],'eps2c')
%%
addpath ../../Functions/

image.all{1} = b_data_das.get_image();
image.tags{1} = 'DAS';
image.all{2} = b_data_mv.get_image();
image.tags{2} = 'MV';
image.all{3} = b_data_cf.get_image();
image.tags{3} = 'CF';
image.all{4} = b_data_pcf.get_image();
image.tags{4} = 'PCF';
image.all{5} = b_data_gcf.get_image();
image.tags{5} = 'GCF';

%%

[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,40,50,scan.x(1)*10^3,scan.x(end)*10^3);

mask=abs(scan.x_axis)<15e-3;
theory=-50*(scan.x_axis(mask)+15e-3)/30e-3;

f89 = figure(89);clf;hold all;
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1});hold on;
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2});
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3});
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4});
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5});
plot(scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');

ylim([-80 0]);
%xlim([-17.5 1.5]);
legend show; grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/gradient'],'eps2c')

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

