clear all;
close all;

%% Download data
url='https://nyhirse.medisin.ntnu.no/ustb/data/gcnr/';   % if not found data will be downloaded from here

%filename='L7_FI_carotid_cross_1.uff';
filename='L7_FI_carotid_cross_2.uff';
tools.download(filename, url, data_path);   

%% Load data
channel_data = uff.channel_data();
channel_data.read([data_path filesep filename],'/channel_data');

channel_data.data = channel_data.data(:,:,:,2)
%% Create Linear Scan 
MLA = 4;
z_axis = linspace(2e-3,35e-3,1500).';
x_axis = linspace(channel_data.sequence(1).source.x,channel_data.sequence(end).source.x,channel_data.N_waves.*MLA); %zeros(channel_data.N_waves.*MLA,1); 
%for n=1:channel_data.N_waves.*MLA
%    x_axis(n)=channel_data.sequence(n).source.x;
    
%end
scan=uff.linear_scan('x_axis',x_axis','z_axis',z_axis);



%% Regions


if strcmp(filename,'L7_FI_carotid_cross_1.uff')
    % cyst geometry -> this should go in the uff
    x0=0.5e-3;                
    z0=14.5e-3; 
    r=2.8e-3;    
    skip=6.5e-3;
    save_path = [ustb_path,filesep,'publications',filesep,'TUFFC',filesep,'Rodriguez-Molares_et_al_Generalized_Contrast_to_Noise_ratio',filesep,'Figures/in_vivo/carotid_cross_1']
elseif strcmp(filename,'L7_FI_carotid_cross_2.uff')

    x0=0.2e-3;                
    z0=15.4e-3; 
    r=3e-3;                 
    skip=-12e-3;
    save_path = [ustb_path,filesep,'publications',filesep,'TUFFC',filesep,'Rodriguez-Molares_et_al_Generalized_Contrast_to_Noise_ratio',filesep,'Figures/in_vivo/carotid_cross_2']
else
    error('Please define appropriate regions');
end
%sca=uff.linear_scan('x_axis', x0 + linspace(-4e-3,9e-3,256).','z_axis', z0 + linspace(-4e-3,4e-3,2.5*256).');

% stand off distance <- based on aperture size
r_off = 0.5e-3;                     % overwrite to handle larger pulse duration

% boundaries
ri=r-r_off;
Ai=pi*ri^2;
d=sqrt((scan.x-x0).^2+(scan.z-z0).^2);
l=sqrt(Ai);


% masks
mask_i=d<ri;
mask_o= ((scan.x>(x0+skip-l/2)).*(scan.x<(x0+skip+l/2)).*(scan.z>(z0-l/2)).*(scan.z<(z0+l/2)))>0;

sum(mask_i)
sum(mask_o)

figure;
subplot(2,1,1)
imagesc(scan.x_axis*1e3, scan.z_axis*1e3, reshape(mask_i,[scan.N_z_axis scan.N_x_axis] )); axis image;
subplot(2,1,2)
imagesc(scan.x_axis*1e3, scan.z_axis*1e3, reshape(mask_o,[scan.N_z_axis scan.N_x_axis] )); axis image;

%% Prepare beamforming
pipe=pipeline();
pipe.scan=scan;
pipe.channel_data = channel_data;

pipe.transmit_apodization.window=uff.window.scanline;
pipe.transmit_apodization.MLA = 4;
pipe.transmit_apodization.MLA_overlap = 2;
pipe.transmit_apodization.f_number = 2;

%pipe.transmit_apodization.minimum_aperture = M*mix.probe.pitch;
%pipe.transmit_apodization.maximum_aperture = M*mix.probe.pitch;

pipe.receive_apodization.window=uff.window.boxcar;
pipe.receive_apodization.f_number = 1;
%pipe.receive_apodization.minimum_aperture = M*mix.probe.pitch;
%pipe.receive_apodization.maximum_aperture = M*mix.probe.pitch;

das=midprocess.das();
das.dimension = dimension.both;
b_das = pipe.go({ das });

%%
f = figure();
b_das.plot(f,'DAS',60); hold on;
tools.plot_circle(x0*1e3,z0*1e3,ri*1e3,'r-');
plot(1e3*(x0+skip+[-l/2 l/2 l/2 -l/2 -l/2]),...
     1e3*(z0+[-l/2 -l/2 l/2 l/2 -l/2]),...
     'g--','Linewidth',2);
saveas(f,[save_path,'_DAS'],'eps2c')
%% DAS

% evaluate contrast
GCNR_das = inVivoGCNR(b_das, mask_o, mask_i, 'DAS')
title(gca,['DAS GCNR = ',num2str(GCNR_das,'%0.3f')]);

%% S-DAS

% compresion
slope=10;
xA=[linspace(0,-25,3) -60+linspace(0,25,3)].';
yA=[linspace(0,-25,3)/slope -60+linspace(0,25,3)/slope].';

f=fit(xA,yA,'smoothingspline')
figure;
xx=linspace(-60,0,100); 
plot(xx,xx,'r--','linewidth',2); hold on; grid on; axis equal;
plot(xx,f(xx)','linewidth',2);
xlim([-60 0]);
ylim([-62 2]);
set(gca,'FontSize', 12);
xlabel('Input [dB]');
ylabel('Output [sB]');
legend('linear','S-curve','Location','NorthWest')

% anechoic
b_sdas = uff.beamformed_data(b_das);
b_sdas.data = 20*log10(abs(b_das.data)./max(abs(b_das.data)));
b_sdas.data = 10.^(reshape(f(b_sdas.data),size(b_sdas.data))/20);



%%
% evaluate contrast
f = figure();
b_sdas.plot(f,'S-DAS',60);
saveas(f,[save_path,'_S_DAS'],'eps2c')
GCNR_sdas = inVivoGCNR(b_sdas, mask_o, mask_i, 'S-DAS')
title(gca,['SDAS GCNR = ',num2str(GCNR_sdas,'%0.3f')]);

%% beamforming on transmit
das.dimension = dimension.transmit;
b_transmit = pipe.go({ das });

%% CF
% beamform
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf.transmit_apodization = das.transmit_apodization;
cf.receive_apodization = das.receive_apodization;
cf.input = b_transmit;
cf_in_vivo = cf.go();

%%
f = figure()
cf_in_vivo.plot(f,'CF',80);
saveas(f,[save_path,'_CF'],'eps2c')
GCNR_CF = inVivoGCNR(cf_in_vivo, mask_o, mask_i, 'CF')
title(gca,['CF GCNR = ',num2str(GCNR_CF)]);

%% PCF

% beamform
pcf = postprocess.phase_coherence_factor();
pcf.center_frequency = 5e6;
pcf.dimension = dimension.receive;
pcf.gamma=1;
pcf.transmit_apodization = das.transmit_apodization;
pcf.receive_apodization = das.receive_apodization;
pcf.input = b_transmit;
pcf_in_vivo = pcf.go();

%%
f = figure();
pcf_in_vivo.plot(f,'PCF',80)
saveas(f,[save_path,'_PCF'],'eps2c')
GCNR_PCF = inVivoGCNR(pcf_in_vivo, mask_o, mask_i, 'PCF')
title(gca,['PCF GCNR = ',num2str(GCNR_PCF)]);

%% GCF
% beamform
gcf = postprocess.generalized_coherence_factor;
gcf.dimension = dimension.receive;
gcf.M0=4;
gcf.transmit_apodization = das.transmit_apodization;
gcf.receive_apodization = das.receive_apodization;
gcf.input = b_transmit;
gcf_in_vivo = gcf.go();
%%
f = figure();
gcf_in_vivo.plot(f,'GCF',80) 
saveas(f,[save_path,'_GCF'],'eps2c')
GCNR_GCF = inVivoGCNR(gcf_in_vivo, mask_o, mask_i, 'GCF')
title(gca,['GCF GCNR = ',num2str(GCNR_GCF)]);

%% DMAS
% process DMAS
dmas=postprocess.delay_multiply_and_sum();
dmas.dimension = dimension.receive;
dmas.transmit_apodization = pipe.transmit_apodization;
dmas.receive_apodization = pipe.receive_apodization;
dmas.input = b_transmit;
dmas.channel_data = channel_data;
dmas_in_vivo = dmas.go();
%%
f = figure();
dmas_in_vivo.plot(f,['DMAS'],60);
saveas(f,[save_path,'_DMAS'],'eps2c')
%%
GCNR_DMAS = inVivoGCNR(dmas_in_vivo, mask_o, mask_i, 'DMAS')
title(gca,['DMASGCNR = ',num2str(GCNR_DMAS)]);

%% SLSC 

das.dimension = dimension.transmit;
das.receive_apodization.window = uff.window.none();
b_transmit_slsc = pipe.go({ das });
% important that we use only M elements, centered around the abscissa of the pixel. 
% Changing that will alter the SNR ratio.

% Set up the SLSC postprocess
slsc = postprocess.short_lag_spatial_coherence();
slsc.receive_apodization = das.receive_apodization;
slsc.dimension = dimension.receive;
slsc.channel_data = channel_data;
slsc.maxM = 14;
slsc.input = b_transmit_slsc;
slsc.K_in_lambda = 1;

Q = slsc.maxM./channel_data.probe.N

%%
slsc_in_vivo = slsc.go()
%Clamping
slsc_in_vivo.data(slsc_in_vivo.data < 0) = 0;

%%
f = figure();
slsc_in_vivo.plot(f,'SLSC');hold on;
saveas(f,[save_path,'_SLSC'],'eps2c')
%caxis([0 0.1])
%tools.plot_circle(x0*1e3,z0*1e3,ri*1e3,'r-');
%plot(1e3*(x0+skip+[-l/2 l/2 l/2 -l/2 -l/2]),...
%     1e3*(z0+[-l/2 -l/2 l/2 l/2 -l/2]),...
%     'g--','Linewidth',2);
%%
GCNR_SLSC = inVivoGCNR(slsc_in_vivo, mask_o, mask_i, 'SLSC')
title(gca,['SLSC GCNR = ',num2str(GCNR_SLSC)]);

%%
f = figure()
subplot(2,1,1);
bar([GCNR_das GCNR_sdas GCNR_CF GCNR_GCF GCNR_PCF GCNR_DMAS GCNR_SLSC])
ylabel('GCNR');
xticks(1:7)
xticklabels({'DAS' 'S-DAS' 'CF' 'GCF' 'PCF' 'DMAS' 'SLSC'})
set(gca,'FontSize',15)
saveas(f,[save_path,'_GCNR'],'eps2c')
