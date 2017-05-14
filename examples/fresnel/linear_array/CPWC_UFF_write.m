% Example of UFF write read routines.

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Last updated: 2017/05/09$

clear all;
close all;

%% PHANTOM
x_sca=[zeros(1,7) -15e-3:5e-3:15e-3];
z_sca=[5e-3:5e-3:35e-3 20e-3*ones(1,7)];
N_sca=length(x_sca);
pha=uff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[x_sca.', zeros(N_sca,1), z_sca.', ones(N_sca,1)];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% PROBE
prb=uff.linear_array();
prb.N=128;                  % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=5000e-6; % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=uff.pulse();
pul.center_frequency=5.2e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
N_plane_waves=15;
angles=linspace(-0.3,0.3,N_plane_waves);
seq=uff.wave();
for n=1:N_plane_waves 
    seq(n)=uff.wave();
    seq(n).probe=prb;
    seq(n).source.azimuth=angles(n);
    seq(n).source.distance=Inf;
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
sca=uff.linear_scan(linspace(-20e-3,20e-3,256).', linspace(0e-3,40e-3,256).');
sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.0;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uff.window.tukey50;
bmf.transmit_apodization.f_number=1.0;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go({process.das_mex() process.coherent_compounding()});
b_data.plot();

% save channel data & beamformed data
uff_file=uff('test02.uff','write');
uff_file.write('/',channel_data,'channel_data');
uff_file.write('/',b_data,'b_data');

clear all;
close all;

% read data from file
uff_file=uff('test02.uff');
uff_file.index([],true);
channel_data=uff_file.read('/channel_data')
b_data=uff_file.read('/b_data')

% show beamformed data
b_data.plot();

%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=b_data.scan;
bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.0;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uff.window.tukey50;
bmf.transmit_apodization.f_number=1.0;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go({process.das_mex() process.coherent_compounding()});
b_data.plot();
