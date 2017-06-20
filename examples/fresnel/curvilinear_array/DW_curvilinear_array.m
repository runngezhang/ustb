% Example of Diverging Waves simulation with the USTB built-in Fresnel simulator

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/23$

clear all;
close all;

%% PHANTOM
pha=uff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[zeros(11,1),  zeros(11,1), linspace(10e-3,160e-3,11).', ones(11,1);...
            linspace(-70e-3,70e-3,11).',  zeros(11,1), 70e-3*ones(11,1), ones(11,1)];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% PROBE
prb=uff.curvilinear_array();
prb.N=128;                  % number of elements 
prb.pitch=508e-6;
prb.element_width=408e-6;
prb.radius=60e-3;
prb.plot(fig_handle);

%% PULSE
pul=uff.pulse();
pul.center_frequency=3.2e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
N=15;                             % number of diverging waves
x0=linspace(-30e-3,30e-3,N);
z0=-prb.radius;
seq=uff.wave();
for n=1:N 
    seq(n)=uff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[x0(n) 0 z0];
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
sca=uff.sector_scan();
sca.azimuth_axis=linspace(-prb.maximum_angle,prb.maximum_angle,256).';
sca.depth_axis=linspace(prb.radius,prb.radius+180e-3,256).';
sca.apex=uff.point('xyz',[0 0 -prb.radius]);
sca.plot(fig_handle,'Scenario');    % show mesh

%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.origo=sca.apex;

% beamforming
b_data=bmf.go({process.das_matlab() process.coherent_compounding()});

% show
h_fig=b_data.plot();
