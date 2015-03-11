%% Computation of a STAI dataset with Field II and beamforming with USTB
%
% This example shows how to load the data from a Field II simulation into a
% STA class, and then demodulate and beamform it with the USTB routines. It
% compared the resulting PSF with the theoterical sinc function.
% The Field II simulation program (field-ii.dk) should be in MATLAB's path.
%
% date:     11.03.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

%% basic constants
c0=1540;     % Speed of sound [m/s]
fs=100e6;   % Sampling frequency [Hz]
dt=1/fs;    % Sampling step [s] 

%% field II initialisation
field_init;
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% transducer definition L9-4/38 Ultrasonix
f0=5e6;             % Transducer center frequency [Hz]
bw=0.1;             % probe bandwidth [1]
lambda=c0/f0;       % Wavelength [m]
height=6e-3;        % Height of element [m]
pitch=0.3048e-3;    % pitch [m]
kerf=0.035e-3;      % gap between elements [m]
width=pitch-kerf;   % Width of element [m]
lens_el=19e-3;      % position of the elevation focus
N_elements=128;     % Number of elements

%% pulse definition
t0=(-1.0/bw/f0): dt : (1.0/bw/f0);
% n_cycles=20;
% n_samples=floor(n_cycles/f0/dt);
% excitation=square(2*pi*f0*(0:(n_samples-1))*dt);
excitation=1;
impulse_response=gauspuls(t0, f0, bw);
two_ways_ir= conv(conv(impulse_response,impulse_response),excitation);
if mod(length(impulse_response),2)
    lag=(length(two_ways_ir)-1)/2;          
else
    lag=(length(two_ways_ir))/2;
end

% show the pulse to check that the lag estimation is on place (and that the pulse is symmetric)
figure;
plot((0:(length(two_ways_ir)-1))*dt -lag*dt,two_ways_ir); hold on; grid on; axis tight
plot((0:(length(two_ways_ir)-1))*dt -lag*dt,abs(hilbert(two_ways_ir)),'r')
plot([0 0],[min(two_ways_ir) max(two_ways_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');

%% aperture objects
% definition of the mesh geometry
noSubAz=round(width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl=round(height/(lambda/8));       % number of subelements in the elevation direction
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
sca=[0 0 20e-3];             % list with the scatterers coordinates [m]
amp=1;                       % list with the scatterers amplitudes
cropat=round(2*40e-3/c0/dt); % maximum time sample, samples after this will be dumped

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
    
    % do calculation
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
sta_dataset=sta('Field II, STA, RF format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c0,...                            % reference speed of sound (m/s)
      t_out.',...                       % time vector (s)
      STA,...                           % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,2)]);      % probe geometry [x, y, z] (m)

% and demodulate  
sta_dataset.demodulate(true);

%% define a reconstruction

% define a scan
scan=linear_scan();
scan.x_axis=linspace(-4e-3,4e-3,256).';               % x vector [m]
scan.z_axis=linspace(16e-3,24e-3,256).';                 % z vector [m]

% define a synthetic orientation
F_number=2;
orien=orientation();
orien.transmit_beam=beam(F_number,E.apodization_type.boxcar);  
orien.receive_beam=beam(F_number,E.apodization_type.boxcar);

% define a reconstruction 
sta_image=reconstruction();
sta_image.scan=scan;
sta_image.orientation=orien;

%% beamform and show
sta_dataset.image_reconstruction(sta_image);
im=sta_image.show();

%% compare lateral profile to sinc
lateral_profile=im(128,:);
lateral_profile=lateral_profile-max(lateral_profile);
theoretical_profile=20*log10(sinc(1/F_number/lambda*scan.x_axis).^2);

figure;
plot(scan.x_axis*1e3,lateral_profile); hold on; grid on; 
plot(scan.x_axis*1e3,theoretical_profile,'r-'); 
legend('Simulation','Theoretical');
xlabel('x [mm]');
ylabel('Amplitude [dB]');

