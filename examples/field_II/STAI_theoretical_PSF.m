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

% Plot the pulse from USTB simulation
pulse.plot([],'2-way pulse for Fresnel simulator');

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

%% PHANTOM
pha=uff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[0,  0, 20e-3, 1];    % point scatterer position [m]
fig_handle=pha.plot();   
cropat=round(1.1*2*sqrt((max(pha.points(:,1))-min(probe.x))^2+max(pha.points(:,3))^2)/c0/dt);   % maximum time sample, samples after this will be dumped

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
    [v,t]=calc_scat_multi(Th, Rh, pha.points(1:3), pha.points(4));
    
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
    
    seq(n).apodization = uff.apodization();
    seq(n).apodization.window=uff.window.sta;
    seq(n).apodization.apex=seq(n).source;
end
close(wb);

%% CHANNEL DATA
channel_data_field_ii = uff.channel_data();
channel_data_field_ii.sampling_frequency = fs;
channel_data_field_ii.sound_speed = c0;
channel_data_field_ii.initial_time = 0;
channel_data_field_ii.pulse = pulse;
channel_data_field_ii.probe = probe;
channel_data_field_ii.sequence = seq;
channel_data_field_ii.data = STA;

%% SCAN
sca=uff.linear_scan(linspace(-4e-3,4e-3,256).', linspace(16e-3,24e-3,256).');
 %% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data_field_ii;
bmf.scan=sca;
bmf.receive_apodization.window=uff.window.boxcar;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uff.window.boxcar;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;
%%
% Delay and sum on receive, then coherent compounding
b_data_field_ii =bmf.go({process.das_mex() process.coherent_compounding()});

%% SIMULATOR
sim=fresnel();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pulse;                  % transmitted pulse
sim.probe=probe;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=channel_data_field_ii.sampling_frequency;  % sampling frequency [Hz]

% we launch the simulation
channel_data_fresnel=sim.go();


%% BEAMFORM data from Fresnel simulation
bmf.channel_data=channel_data_fresnel;
% Delay and sum on receive, then coherent compounding
b_data_fresnel =bmf.go({process.das_mex() process.coherent_compounding()});


%% Display images
figure(101);
ax1 = subplot(121);
ax2 = subplot(122);

b_data_field_ii.plot(ax1,'Field II Simulation')
b_data_fresnel.plot(ax2,'Fresnel Simulation')


%% compare lateral profile to sinc
img_field_ii = b_data_field_ii.get_image;
lateral_profile_field_ii=img_field_ii(128,:);
lateral_profile_field_ii=lateral_profile_field_ii-max(lateral_profile_field_ii);

img_fresnel = b_data_fresnel.get_image;
lateral_profile_fresnel=img_fresnel(128,:);
lateral_profile_fresnel=lateral_profile_fresnel-max(lateral_profile_fresnel);

theoretical_profile=20*log10(sinc(1/bmf.receive_apodization.f_number(1)/lambda*b_data_field_ii.scan.x_axis).^2);

figure;
plot(b_data_field_ii.scan.x_axis*1e3,lateral_profile_field_ii); hold all; grid on; 
plot(b_data_field_ii.scan.x_axis*1e3,lateral_profile_fresnel,'k'); 
plot(b_data_field_ii.scan.x_axis*1e3,theoretical_profile,'r'); 
legend('Field II Simulation','Fresnel Simulation','Theoretical');
xlabel('x [mm]');
ylabel('Amplitude [dB]');
title('Lateral (x-axis) profile ');

%% compare axial profile
axial_profile_field_ii=img_field_ii(:,128);
axial_profile_field_ii=axial_profile_field_ii-max(axial_profile_field_ii);

axial_profile_fresnel=img_fresnel(:,128);
axial_profile_fresnel=axial_profile_fresnel-max(axial_profile_fresnel);

figure;
plot(b_data_field_ii.scan.x_axis*1e3,axial_profile_field_ii); hold all; grid on; 
plot(b_data_field_ii.scan.x_axis*1e3,axial_profile_fresnel,'k'); 
legend('Field II Simulation','Fresnel Simulation');
xlabel('z [mm]');
ylabel('Amplitude [dB]');
title('Axial (z-axis) profile ');

