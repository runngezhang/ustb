% Example of Conventional Focused Imaging (single focal depth)
% for a phased-array and sector scan. Simulated with the USTB's Fresnel
% simulator 
%
% This example uses a dummy example with multiple frames just moving a
% single point scatterer 1 mm between frames.

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal (olemarius@olemarius.net)
%   $Date: 2017/03/11$

clear all;
close all;

%% PHANTOM
N_sca=1;                           % number of scatterers
x_sca=-20e-3;
z_sca=30e-3;
p=[x_sca zeros(N_sca,1) z_sca];
N_frames=5;                         % number of frames
N_beams=128;                        % number of focused beams
alpha=35*pi/180;                    % scatterer direction [rad]
v_mag=0.25;                         % scatterer velocity magnitude [m/s]
v=v_mag*ones(N_sca,1)*[cos(alpha) 0 sin(alpha)]; % scatterer velocity [m/s m/s m/s]
PRF=4000;                           % pulse repetition frequency [Hz]
fig_handle=figure();
for n=1:N_frames*N_beams
    pha(n)=uff.phantom();
    pha(n).sound_speed=1540;            % speed of sound [m/s]
    pha(n).points=[p+v*(n-1)/PRF, ones(N_sca,1)];    % point scatterer position [m]
    pha(n).plot(fig_handle);             
end

%% PROBE
prb=uff.linear_array();
prb.N=64;                   % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=7000e-6; % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=uff.pulse();
pul.center_frequency=3e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;   % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
azimuth_axis=linspace(-35*pi/180,35*pi/180,N_beams).';
depth=40e-3;
seq=uff.wave();
for n=1:N_beams
    seq(n)=uff.wave();
    seq(n).probe=prb;
    
    seq(n).source=uff.point();
    seq(n).source.azimuth=azimuth_axis(n);
    seq(n).source.distance=depth;
    
    seq(n).apodization.window=uff.window.tukey50;
    seq(n).apodization.f_number=1.7;
    seq(n).apodization.scan.xyz=seq(n).source.xyz;
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% SIMULATOR
sim=fresnel();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% SCAN
depth_axis=linspace(5e-3,80e-3,256).';
sca=uff.sector_scan();
for n=1:N_beams
    sca(n)=uff.sector_scan(azimuth_axis(n),depth_axis);
    sca(n).plot(fig_handle,'Scenario');    
end
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.7;

% beamforming
b_data=bmf.go({process.das_mex() process.stack()});


%% show
b_data.plot();