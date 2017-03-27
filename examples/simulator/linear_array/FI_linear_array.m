% Example of Conventional Focused Imaging (single focal depth)
% for a linear-array and linear scan. Simulated with the USTB's Fresnel
% simulator 

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
N=200;                      % number of focused beams
x_axis=linspace(-2e-3,2e-3,N).';
z0=40e-3;
seq=huff.wave();
for n=1:N 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[x_axis(n) 0 z0];
    
    seq(n).apodization.window=huff.window.tukey50;
    seq(n).apodization.f_number=1.7;
    seq(n).apodization.apex.distance=Inf;
    seq(n).apodization.scan.xyz=seq(n).source.xyz;
    
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
 
%% SCAN
z_axis=linspace(39e-3,41e-3,100).';
sca=huff.linear_scan();
for n=1:length(x_axis)
    sca(n)=huff.linear_scan(x_axis(n),z_axis);
    sca(n).plot(fig_handle,'Scenario');    
end
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
bmf.receive_apodization.window=huff.window.tukey50;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go(@bmf.matlab,@postprocess.stack);

% show
b_data.plot();