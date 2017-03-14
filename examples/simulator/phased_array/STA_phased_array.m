% Example of Synthetic Transmit Aperture simulation with USTB's Fresnel simulator
% for a phased array and a SECTOR_SCAN

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/14$

clear all;
close all;

%% PHANTOM
pha=huff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[0,  0, 40e-3, 1];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% PROBE
prb=huff.linear_array();
prb.N=64;                   % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=7000e-6; % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=huff.pulse();
pul.center_frequency=3e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;   % fractional bandwidth [unitless]
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
sim.phantom=pha;                                % phantom
sim.pulse=pul;                                  % transmitted pulse
sim.probe=prb;                                  % probe
sim.sequence=seq;                               % beam sequence
sim.sampling_frequency=4*pul.center_frequency;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% SCAN
sca=huff.sector_scan(linspace(-10*pi/180,10*pi/180,200).', linspace(35e-3,45e-3,100).');
sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
bmf.receive_apodization.window=huff.window.tukey50;
bmf.receive_apodization.f_number=1.7;
bmf.transmit_apodization.window=huff.window.tukey50;
bmf.transmit_apodization.f_number=1.7;

% beamforming
b_data=bmf.go(@postprocess.coherent_compound);

% show
b_data.plot();

