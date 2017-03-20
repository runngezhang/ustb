% Example of movement simulation with USTB's Fresnel simulator
% using Coherent Plane-Wave Compounding sequence on a linear array.

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/13$

clear all;
close all;

%% PHANTOM
alpha=-45*pi/180;                 % vessel direction [rad]
N_sca=1;                          % number of scatterers
%x_sca=random('unif',-10e-3,10e-3,N_sca,1);
%z_sca=random('unif',15e-3,25e-3,N_sca,1);
x_sca=-1e-3;
z_sca=21e-3;
p=[x_sca zeros(N_sca,1) z_sca+x_sca*sin(alpha)];
v=0.9754*ones(N_sca,1)*[cos(alpha) 0 sin(alpha)]; % scatterer velocity [m/s m/s m/s]
PRF=10000;                           % pulse repetition frequency [Hz]
N_plane_waves=3;                     % number of plane wave
N_frames=10;                         % number of frames
fig_handle=figure();
for n=1:N_plane_waves*N_frames
    pha(n)=huff.phantom();
    pha(n).sound_speed=1540;            % speed of sound [m/s]
    pha(n).points=[p+v*(n-1)/PRF, ones(N_sca,1)];    % point scatterer position [m]
    pha(n).plot(fig_handle);             
end
             
%% PROBE
prb=huff.linear_array();
prb.N=128;                  % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=5000e-6; % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=huff.pulse();
pul.center_frequency=5.2e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
angles=linspace(-0.3,0.3,N_plane_waves);
seq=huff.wave();
for n=1:N_plane_waves 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.azimuth=angles(n);
    seq(n).source.distance=Inf;
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.PRF=PRF;                    % pulse repetition frequency [Hz]
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% SCAN
sca=huff.linear_scan(linspace(-19e-3,19e-3,256).', linspace(0e-3,40e-3,256).');
sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
bmf.receive_apodization.window=huff.window.tukey50;
bmf.receive_apodization.f_number=1.0;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=huff.window.tukey50;
bmf.transmit_apodization.f_number=1.0;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go(@bmf.matlab,@postprocess.coherent_compound);

%% show
b_data.plot([],['Beamformed data'],40);