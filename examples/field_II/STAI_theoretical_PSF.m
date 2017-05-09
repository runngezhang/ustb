%% Computation of a STAI dataset with Field II and beamforming with USTB
%
% This example shows how to load the data from a Field II simulation into 
% USTB objects, and then beamform it with the USTB routines. It compared 
% the resulting PSF with the theoterical sinc function.
% The Field II simulation program (field-ii.dk) should be in MATLAB's path.
%
% date:     11.03.2015
% updated:  09.05.2017
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%           Ole Marius Hoel Rindal <olemarius@olemarius.net>

clear all;
close all;

%% basic constants
c0=1540;     % Speed of sound [m/s]
fs=100e6;    % Sampling frequency [Hz]
dt=1/fs;     % Sampling step [s] 

%% field II initialisation
field_init;
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% transducer definition L9-4/38 Ultrasonix
probe = uff.linear_array();
f0                      =5e6;             % Transducer center frequency [Hz]
lambda                  =c0/f0;           % Wavelength [m]
probe.element_height    =6e-3;            % Height of element [m]
probe.pitch             =0.3048e-3;       % probe.pitch [m]
kerf                    =0.035e-3;        % gap between elements [m]
probe.element_width     =probe.pitch-kerf;% Width of element [m]
lens_el                 =19e-3;           % position of the elevation focus
probe.N                 =128;             % Number of elements

%% pulse definition
pulse = uff.pulse(f0);
pulse.fractional_bandwidth = 0.1;             % probe bandwidth [1]
t0=(-1.0/pulse.fractional_bandwidth /f0): dt : (1.0/pulse.fractional_bandwidth /f0);
excitation=1;
impulse_response=gauspuls(t0, f0, pulse.fractional_bandwidth );
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
noSubAz=round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl=round(probe.element_height/(lambda/8));       % number of subelements in the elevation direction
Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 

% setting excitation, impulse response and baffle
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);

%% phantom
sca=[0 0 20e-3];             % list with the scatterers coordinates [m]
amp=1;                       % list with the scatterers amplitudes
cropat=round(2*40e-3/c0/dt); % maximum time sample, samples after this will be dumped

%% output data
t_out=0:dt:((cropat-1)*dt);                 % output time vector
STA=zeros(cropat,probe.N,probe.N);    % impulse response channel data
%% Compute STA signals
disp('Field II: Computing STA dataset');
wb = waitbar(0, 'Field II: Computing STA dataset');
for n=1:probe.N
    waitbar(n/probe.N, wb);

    % transmit aperture
    xdc_apodization(Th, 0, [zeros(1,n-1) 1 zeros(1,probe.N-n)]);
    xdc_focus_times(Th, 0, zeros(1,probe.N));
    
    % receive aperture    
    xdc_apodization(Rh, 0, ones(1,probe.N));
    xdc_focus_times(Rh, 0, zeros(1,probe.N));
    
    % do calculation
    [v,t]=calc_scat_multi(Th, Rh, sca, amp);
    
    % lag compensation
    t_in=(0:dt:((size(v,1)-1)*dt))+t-lag*dt + probe.r(n)/c0;
    v_aux=interp1(t_in,v,t_out,'linear',0);

    % build the dataset
    STA(:,:,n)=v_aux;
    
    %% SEQUENCE GENERATION
    seq(n)=uff.wave();
    seq(n).probe=probe;
    seq(n).source.xyz=[probe.x(n) probe.y(n) probe.z(n)];
    seq(n).sound_speed=c0;
end
close(wb);

%% CHANNEL DATA
channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.data = STA;

%% SCAN
sca=uff.linear_scan(linspace(-4e-3,4e-3,256).', linspace(16e-3,24e-3,256).');

%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
bmf.receive_apodization.window=uff.window.boxcar;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uff.window.boxcar;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;

% Delay and sum on receive, then coherent compounding
b_data=bmf.go({process.das_mex() process.coherent_compounding()});
% Display image
b_data.plot()

%% compare lateral profile to sinc
im = b_data.get_image;
lateral_profile=im(128,:);
lateral_profile=lateral_profile-max(lateral_profile);
theoretical_profile=20*log10(sinc(1/bmf.receive_apodization.f_number(1)/lambda*b_data.scan.x_axis).^2);

figure;
plot(b_data.scan.x_axis*1e3,lateral_profile); hold on; grid on; 
plot(b_data.scan.x_axis*1e3,theoretical_profile,'r-'); 
legend('Simulation','Theoretical');
xlabel('x [mm]');
ylabel('Amplitude [dB]');

