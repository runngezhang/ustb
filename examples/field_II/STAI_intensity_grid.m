%% Measuring the output dynamic range of adaptive beamforming algorithms
%
% date:     11.03.2016
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

clear all;
close all;

%% basic constants
c0=1540;     % Speed of sound [m/s]
fs=100e6;   % Sampling frequency [Hz]
dt=1/fs;    % Sampling step [s] 

%% field II initialisation
field_init;
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% transducer definition P4-2v Verasonics
f0=2.5e6;           % Transducer center frequency [Hz]
bw=0.67;            % probe bandwidth [1]
lambda=c0/f0;       % Wavelength [m]
height=5e-3;        % Height of element [m]
pitch=0.300e-3;     % pitch [m]
kerf=0.050e-3;      % gap between elements [m]
width=pitch-kerf;   % Width of element [m]
lens_el=60e-3;      % position of the elevation focus
N_elements=64;     % Number of elements
pulse_duration=2.5; % pulse duration [cycles]

%% pulse definition
t0 = (-1/bw/f0): dt : (1/bw/f0);
impulse_response = gauspuls(t0, f0, bw);
te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2;   

% show the pulse to check that the lag estimation is on place (and that the pulse is symmetric)
figure;
plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');

%% aperture objects
% definition of the mesh geometry
noSubAz=round(width/(lambda/8));              % number of subelements in the azimuth direction
noSubEl=round(height/(lambda/8));              % number of subelements in the elevation direction
Th = xdc_linear_array (N_elements, width, height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
Rh = xdc_linear_array (N_elements, width, height, kerf, noSubAz, noSubEl, [0 0 Inf]); 

% setting excitation, impulse response and baffle
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);

%% get geometrical center of elements
data = xdc_get(Th,'rect');
geo=data(24:26,:);
x0=zeros(1,N_elements);
for n=1:N_elements
    n_ini=(n-1)*noSubAz*noSubEl+1;
    n_fin=n_ini+noSubAz*noSubEl-1;
    x0(n)=mean(geo(1,n_ini:n_fin));
end

%% phantom
x_point=-60e-3:10e-3:60e-3;                                 % x-coordinate of the scatterers [m]
z_point=10e-3:10e-3:90e-3;                                  % z-coordinate of the scatterers [m]
[xxp zzp]=meshgrid(x_point,z_point);                        
mask=abs(atan2(xxp,zzp-10e-3))<=pi/4;
mask=mask&(sqrt(xxp.^2+(zzp-10e-3).^2)<90e-3);
xxp=xxp(mask);
zzp=zzp(mask);
N_sca=length(zzp(:));                                       % total number of scatterers
sca=[xxp(:) zeros(N_sca,1) zzp(:)];                         % list with the scatterers coordinates [m]

figure;
plot(sca(:,1),sca(:,3),'b.'); hold on; axis equal; grid on;

% compensate for energy divergency 
att = inline('2307*x.^2 -570*x+ 5.803');
amp_dB=round(random('unif',-60,0,N_sca,1))-att(sqrt(sum(sca.^2,2)));              % list with the scatterers amplitudes

amp=sqrt(10.^(amp_dB/10));
cropat=round(1.1*2*sqrt((max(sca(:,1))-min(x0))^2+max(sca(:,3))^2)/c0/dt);   % maximum time sample, samples after this will be dumped

%% output data
t_out=0:dt:((cropat-1)*dt);                 % output time vector
STA=zeros(cropat,N_elements,N_elements);    % impulse response channel data

%% Compute STA signals
disp('Field II: Computing STA dataset');
wb = waitbar(0, 'Field II: Computing STA dataset');
for n=1:N_elements
    waitbar(n/N_elements, wb);

    % transmit aperture
    xdc_apodization(Th, 0, [zeros(1,n-1) 1 zeros(1,N_elements-n)]);
    xdc_focus_times(Th, 0, zeros(1,N_elements));
    
    % receive aperture    
    xdc_apodization(Rh, 0, ones(1,N_elements));
    xdc_focus_times(Rh, 0, zeros(1,N_elements));
    
    % do calculation -> include TGC in the simulation
    [v,t]=calc_scat_multi(Th, Rh, sca, amp);
    
    % lag compensation
    t_in=(0:dt:((size(v,1)-1)*dt))+t-lag*dt;
    v_aux=interp1(t_in,v,t_out,'linear',0);

    % build the dataset
    STA(:,:,n)=v_aux;
end
close(wb);

%% USTB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define a sta object
dataset=sta('Field II, STA, RF format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c0,...                            % reference speed of sound (m/s)
      t_out.',...                       % time vector (s)
      STA,...                           % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,2)]);      % probe geometry [x, y, z] (m)

% and demodulate  
dataset.demodulate(true);

% reconstruction
recons=reconstruction();

% scan
angle_span=90*pi/180;
recons.scan=sector_scan();
recons.scan.depth_axis=linspace(0,100e-3,512).';
recons.scan.azimuth_axis=linspace(-angle_span/2,angle_span/2,512).';

% orientation
recons.orientation=orientation();
recons.orientation.transmit_beam=beam(0,E.apodization_type.none);
recons.orientation.receive_beam=beam(0,E.apodization_type.none);

% beamform and show
dataset.image_reconstruction(recons);
recons.show();

%% delay
w0=2*pi*dataset.modulation_frequency;
sig=zeros(N_sca,dataset.channels,dataset.firings,dataset.frames);
wb=waitbar(0,'Delaying');
distance=sqrt(xxp.^2+zzp.^2);
for f=1:dataset.frames
    for ntx=1:dataset.channels
        waitbar(((f-1)*dataset.channels+ntx)/dataset.channels/dataset.frames);
        TF=sqrt((dataset.geom(ntx,1)-xxp).^2+(dataset.geom(ntx,3)-zzp).^2);
        for nrx=1:dataset.channels
            RF=sqrt((dataset.geom(nrx,1)-xxp).^2+(dataset.geom(nrx,3)-zzp).^2);
            delay=(RF+TF)/dataset.c0;
            phase_shift=exp(1i.*w0*delay-1i.*w0*distance/dataset.c0); 
            sig(:,nrx,ntx,f)=sig(:,nrx,ntx,f)+phase_shift.*interp1(dataset.time,dataset.data(:,nrx,ntx,f),delay,'linear',0);
        end
    end
end
close(wb);

%% and sum
das=sum(sig(:,:),2);
DAS=abs(das);DAS=20*log10(DAS./max(DAS(:)));

%% Mallart-Fink coherence factor 
% R. Mallart and M. Fink, “Adaptive focusing in scattering media through sound-speed inhomogeneities: 
% The van Cittert Zernike approach and focusing criterion,” J. Acoust. Soc. Am., vol. 96, no. 6, pp. 3721–3732, 1994
coherent=abs(sum(sig(:,:),2)).^2;
incoherent=sum(abs(sig(:,:).^2),2);
cf=coherent./dataset.channels^2./incoherent;
CF=abs(cf.*das); CF=20*log10(CF./max(CF(:)));

%% Camacho-Fritsch Absolute Phase Coherence factor
% http://cdn.intechopen.com/pdfs-wm/14874.pdf
sigma_0=pi/sqrt(3);
apca=1-std(angle(sig(:,:)),[],2)/sigma_0;
apca(apca<0)=0;
APCA=abs(apca.*das); APCA=20*log10(APCA./max(APCA(:)));

%% show results
figure;
plot(DAS,DAS,'k-'); axis equal; hold on; grid on;
plot(DAS,CF,'r+'); 
plot(DAS,APCA,'bo'); 
xlabel('Input dynamic range');
ylabel('Output dynamic range');

legend('DAS','Mallart-Fink CF', 'Camacho-Fritsch APCA');





