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


%% SLSC using M elements (hacking the SLSC postprocess)

% important that we use only M elements, centered around the abscissa of the pixel. 
% Changing that will alter the SNR ratio.
das.dimension = dimension.transmit;
b_data_tx = pipe.go({das});

%%
% Set up the SLSC postprocess
slsc = postprocess.short_lag_spatial_coherence();
slsc.receive_apodization = das.receive_apodization;
slsc.dimension = dimension.receive;
slsc.channel_data = mix;
slsc.maxM = 14;
slsc.input = b_data_tx;
slsc.K_in_lambda = 1;

Q = slsc.maxM./M

% Pick out the M channels that contains data, thus only the M elemetns
% centered around the abscissa of the pixel. This is taken care of by the
% apodization so we can just pick up the element signals that are larger than 0.

% Do you want to normalize the SLSC values??
normalize_slsc = 1;

% Data buffers
aux_data = zeros(b_data_tx.scan.N_z_axis*b_data_tx.scan.N_x_axis,1,1,mix.N_frames);
aux_data_clamped = zeros(b_data_tx.scan.N_z_axis*b_data_tx.scan.N_x_axis,1,1,mix.N_frames);
data_cube_M_elements = complex(zeros(b_data_tx.scan.N_z_axis,b_data_tx.scan.N_x_axis,M,1));
for f = 1:mix.N_frames
    f
    % Reshape the beamformed data as a cube (Z,X,Elements)
    data_cube = reshape(b_data_tx.data(:,:,1,f),sca.N_z_axis,sca.N_x_axis,mix.N_channels);
    for x = 1: b_data_tx.scan.N_x_axis
        sum_over_z = abs(sum(squeeze(data_cube(:,x,:)),1));
        elements_with_data = sum_over_z>0;
        data_cube_M_elements(:,x,:) = data_cube(:,x,elements_with_data);
    end
    
    % Call the actual implementation of the SLSC calculations, but using
    % the cube which has only M elements
    [image,slsc_values] = slsc.short_lag_spatial_coherence_implementation(data_cube_M_elements);
    image(image<0) = 0; % Set negative coherence values to zero
    
    slsc_img = squeeze(sum(slsc_values(:,:,:),2));
    
    % Make one clamped version
    slsc_img_clamped = slsc_img;
    slsc_img_clamped(slsc_img_clamped < 0 ) = 0;
    aux_data_clamped(:,1,1,f) = slsc_img_clamped(:);
    
    % In the other we  can shift the negative coherence to something
    % positive.
    slsc_img = slsc_img + abs(min(slsc_img(:)));
    aux_data(:,1,1,f) = slsc_img(:);
end

% Put the resulting SLSC images in a beamformed data
b_slsc_M = uff.beamformed_data();
b_slsc_M.scan = sca;
b_slsc_M.data = aux_data;

b_slsc_M_clamped = uff.beamformed_data();
b_slsc_M_clamped.scan = sca;
b_slsc_M_clamped.data = aux_data_clamped;

%%
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_slsc_M_clamped, mask_o, mask_i, 'SLSC M = 55 clamped');
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_slsc_M, mask_o, mask_i, 'SLSC M = 55 shifted');
%% Plotting all the SLSC values in one plot for comparison
figure(123);
subplot(1,2,1+normalize_slsc)
plot(10*log10(channel_SNR),GCNR_N1,'MarkerSize',7,'DisplayName','N = 1'); hold on;
plot(10*log10(channel_SNR),GCNR_N5,'MarkerSize',7,'DisplayName','N = 5'); hold on;
plot(10*log10(channel_SNR),GCNR_N10,'MarkerSize',7,'DisplayName','N = 10'); hold on;
plot(10*log10(channel_SNR),GCNR_N20,'MarkerSize',7,'DisplayName','N = 20'); hold on;
plot(10*log10(channel_SNR),GCNR_N30,'MarkerSize',7,'DisplayName','N = 30'); hold on;
plot(10*log10(nunu),GCNR0(C0(nunu)),'r--','linewidth',2,'DisplayName','Eq.(34)'); hold on; grid on; axis tight square;
set(gca,'FontSize', 12);
ylim([0 1])
xlabel('10 log_{10} \nu_S/\nu_N');
ylabel('GCNR');
legend('show','Location','SouthEast')
title(['SLSC M=55 ','normalized = ',num2str(normalize_slsc)]);
set(gca,'FontSize', 14);



%% SLSC with full aperture
das.dimension = dimension.transmit;
pipe.transmit_apodization.window = uff.window.none;
pipe.receive_apodization.window = uff.window.none;

% Set up the SLSC postprocess
slsc = postprocess.short_lag_spatial_coherence();
slsc.receive_apodization = pipe.receive_apodization;
slsc.dimension = dimension.receive;
slsc.channel_data = mix;
slsc.maxM = 20;
slsc.K_in_lambda = 1;

b_slsc = pipe.go({das slsc});

%% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_slsc, mask_o, mask_i, 'SLSC (N=20) full aperture');


