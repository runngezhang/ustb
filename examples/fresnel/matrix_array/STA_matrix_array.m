% Example of Synthetic Transmit Aperture simulation with USTB built-in Fresnel simulator
% using a 16 x 16 2D matrix probe

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/25$

clear all;
close all;

F_number=1.0;
z0=10e-3;

%% PHANTOM
pha=uff.phantom();
pha.sound_speed=1540;                 % speed of sound [m/s]
pha.points=[0e-3, 0, z0-2e-3, 1;...        % point scatterer positions [m]
            0, -2e-3, z0, 1;...    
            -3e-3, 0, z0+2e-3, 1];      
fig_handle=pha.plot();             
             
%% PROBE
prb=uff.matrix_array();
prb.N_x=32;                 % number of elements 
prb.N_y=4;                  % number of elements 
prb.pitch_x=400e-6;         % probe pitch in azimuth [m]
prb.pitch_y=400e-6;         % probe pitch in azimuth [m]
prb.element_width=570e-6;   % element width [m]
prb.element_height=570e-6;  % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=uff.pulse();
pul.center_frequency=3e6;         % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
seq=uff.wave();
for n=1:prb.N_elements 
    seq(n)=uff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[prb.x(n) prb.y(n) prb.z(n)];
    
    seq(n).apodization.window=uff.window.sta;
    seq(n).apodization.apex=seq(n).source;
    
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
sca=uff.linear_3D_scan(linspace(-6.4e-3,6.4e-3,200).', linspace(z0-3.2e-3,z0+3.2e-3,100).',0);
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
b_data=bmf.go({process.das_matlab process.coherent_compounding});

% show
fig_plot=pha.plot();
b_data.plot(fig_plot); 


