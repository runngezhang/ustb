clear all;
close all;

%% Download data
url='https://nyhirse.medisin.ntnu.no/ustb/data/gcnr/';   % if not found data will be downloaded from here
filename='insilico_20.uff';
tools.download(filename, url, data_path);   

%% Load data
mix = uff.channel_data();
mix.read([data_path filesep filename],'/mix');
channel_SNR = h5read([data_path filesep filename],'/channel_SNR');;

%% Scan
sca=uff.linear_scan('x_axis',linspace(-6e-3,6e-3,256).','z_axis', linspace(14e-3,26e-3,2.5*256).');

%% Regions

% cyst geometry -> this should go in the uff
x0=0e-3;                
z0=20e-3; 
r=3e-3;                 

% stand off distance <- based on aperture size
M = 55;                             % aperture size
aperture = M * mix.probe.pitch;     % active aperture
F = z0 / aperture;                  % F-number
r_off = round(1.2 * mix.lambda * F, 5); % stand-off distance (rounded to 0.01 mm) 

% boundaries
ri=r-r_off;
ro=r+r_off;
rO=sqrt(ri^2+ro^2);
Ai=pi*ri^2;
Ao=pi*rO^2-pi*ro^2;
d=sqrt((sca.x-x0).^2+(sca.z-z0).^2);

% masks
mask_i=d<ri;
mask_o=(d>ro)&(d<rO);

%% Prepare beamforming
pipe=pipeline();
pipe.scan=sca;
pipe.channel_data = mix;

pipe.transmit_apodization.window=uff.window.flat;
pipe.transmit_apodization.f_number = 1;
pipe.transmit_apodization.minimum_aperture = M*mix.probe.pitch;
pipe.transmit_apodization.maximum_aperture = M*mix.probe.pitch;

pipe.receive_apodization.window=uff.window.flat;
pipe.receive_apodization.f_number = 1;
pipe.receive_apodization.minimum_aperture = M*mix.probe.pitch;
pipe.receive_apodization.maximum_aperture = M*mix.probe.pitch;

das=midprocess.das();

%% DAS

% beamform
das.dimension = dimension.both;
b_das = pipe.go({ das });

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_das, mask_o, mask_i, 'DAS');

%% S-DAS

% beamform
das.dimension = dimension.both;

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

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_sdas, mask_o, mask_i, 'S-DAS');

%% CF

% beamform
das.dimension = dimension.transmit;
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf_anechoic = pipe.go({ das cf });

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, cf.CF, mask_o, mask_i, 'CF');

%% PCF

% beamform
das.dimension = dimension.transmit;
pcf = postprocess.phase_coherence_factor();
pcf.center_frequency = 5e6;
pcf.dimension = dimension.receive;
pcf.gamma=1;
pipe.channel_data=mix;
pcf_anechoic = pipe.go({ das pcf });

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, pcf.FCC, mask_o, mask_i, 'PCF');

%% GCF

% beamform
das.dimension = dimension.transmit;
gcf = postprocess.generalized_coherence_factor;
gcf.dimension = dimension.receive;
gcf.M0=4;
pipe.channel_data=mix;
gcf_anechoic = pipe.go({ das gcf });

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, gcf.GCF, mask_o, mask_i, 'GCF');

%% DMAS

% beamform
das.dimension = dimension.transmit;
b_data_tx = pipe.go({das});

% process DMAS
dmas=postprocess.delay_multiply_and_sum();
dmas.dimension = dimension.receive;
dmas.transmit_apodization = pipe.transmit_apodization;
dmas.receive_apodization = pipe.receive_apodization;
dmas.input = b_data_tx;
dmas.channel_data = mix;
b_dmas = dmas.go();

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_dmas, mask_o, mask_i, 'DMAS');


%% SLSC

% important that we use only M elements, centered around the abscissa of the pixel. 
% Changing that will alter the SNR ratio.
das.dimension = dimension.transmit;
pipe.transmit_apodization.window = uff.window.none;
pipe.receive_apodization.window = uff.window.none;
b_data_tx = pipe.go({das});

%%

slsc = postprocess.short_lag_spatial_coherence();
slsc.receive_apodization = das.receive_apodization;
slsc.dimension = dimension.receive;
slsc.channel_data = mix;
slsc.maxM = 10;
slsc.input = b_data_tx;
slsc.K_in_lambda = 1;
slsc_data = slsc.go();
%%
% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, slsc_data, mask_o, mask_i, 'SLSC');

%%
slsc_data.plot([],['SLSC'],[0 1],'none');

