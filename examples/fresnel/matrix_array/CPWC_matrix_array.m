% Example of Coherent Plane-Wave Imaging with USTB built-in Fresnel simulator
% using a 2D matrix probe

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/29$

clear all;
close all;

F_number=1.7;

%% PHANTOM
pha=uff.phantom();
pha.sound_speed=1540;               % speed of sound [m/s]
pha.points=[0, 0, 20e-3, 1];...     % point scatterer position [m]
fig_handle=pha.plot();             
             
%% PROBE
prb=uff.matrix_array();
prb.N_x=16;                 % number of elements 
prb.N_y=16;                 % number of elements 
prb.pitch_x=600e-6;         % probe pitch in azimuth [m]
prb.pitch_y=600e-6;         % probe pitch in azimuth [m]
prb.element_width=570e-6;   % element width [m]
prb.element_height=570e-6;  % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=uff.pulse();
pul.center_frequency=3e6;         % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
alpha_max=1/2/F_number;
N=15;
tx_angles=linspace(-alpha_max,alpha_max,N);
seq=uff.wave();
for n=1:N
    seq(n)=uff.wave();
    seq(n).probe=prb;
    
    seq(n).source.azimuth=tx_angles(n);
    seq(n).source.distance=Inf;
    
    seq(n).apodization.window=uff.window.none;
    
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
sca=uff.linear_3D_scan(linspace(-4e-3,4e-3,200).', linspace(18e-3,22e-3,100).',0);
sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=F_number;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uff.window.tukey50;
bmf.transmit_apodization.f_number=F_number;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go({ process.das_matlab, process.coherent_compounding});

% show
b_data.plot();

