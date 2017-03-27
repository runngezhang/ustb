% Example of Synthetic Transmit Aperture simulation with USTB built-in Fresnel simulator

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/11$

clear all;
%close all;

%% PHANTOM
pha=huff.phantom();
pha.sound_speed=1540;               % speed of sound [m/s]
pha.points=[0, -2e-3, 20e-3, 1];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% PROBE
prb=huff.matrix_array();
prb.N_x=16;                 % number of elements 
prb.N_y=16;                 % number of elements 
prb.pitch_x=600e-6;         % probe pitch in azimuth [m]
prb.pitch_y=600e-6;         % probe pitch in azimuth [m]
prb.element_width=570e-6;   % element width [m]
prb.element_height=570e-6;  % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=huff.pulse();
pul.center_frequency=3e6;         % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
seq=huff.wave();
for n=1:prb.N_elements 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[prb.x(n) prb.y(n) prb.z(n)];
    
    seq(n).apodization.window=huff.window.sta;
    seq(n).apodization.apex=seq(n).source;
    
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
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();

% %% check channel data
% n_beam=1;
% channel_data.plot(n_beam); hold on;
% dtx=norm([prb.x(n_beam)-pha.points(1) prb.y(n_beam)-pha.points(2) prb.z(n_beam)-pha.points(3)])+norm([prb.x(n_beam) prb.y(n_beam) prb.z(n_beam)]);
% drx=sqrt(sum(([prb.x-pha.points(1) prb.y-pha.points(2) prb.z-pha.points(3)]).^2,2));
% delay=(dtx+drx)/pha.sound_speed;
% plot(1:prb.N_elements,delay*1e6,'r--','linewidth',2)
 
%% SCAN
sca=huff.linear_3D_scan(linspace(-4e-3,4e-3,200).', linspace(18e-3,22e-3,100).',pi/2);
sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=huff.window.tukey50;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=huff.window.tukey50;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go(@bmf.matlab,@postprocess.coherent_compound);

% show
b_data.plot();

