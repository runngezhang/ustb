%% Adaptive beamforming in a CPWC simulation
%
% In this example we show how to use adaptive beamforming with the USTB
%
% Author: Ole Marius Hoel Rindal <olemarius@olemarius.net> 
%
% Date : 2017/04/07

clear all;
close all;

%% Phantom

pha=uff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[0,  0, 40e-3, 1];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% Probe

prb=uff.linear_array();
prb.N=128;                  % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=5000e-6; % element height [m]
prb.plot(fig_handle);

%% Pulse
pul=uff.pulse();
pul.center_frequency=5.2e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% Sequence generation

N=3;                           % number of plane waves
angles=linspace(-0.3,0.3,N)   % angle vector [rad]
seq=uff.wave();
for n=1:N 
    seq(n)=uff.wave();
    
    seq(n).source.azimuth=angles(n);
    seq(n).source.distance=Inf;
    
    seq(n).probe=prb;
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% The Fresnel simulator
sim=fresnel();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% Scan
sca=uff.linear_scan(linspace(-2e-3,2e-3,200).', linspace(39e-3,41e-3,100).');
sca.plot(fig_handle,'Scenario');    % show mesh
 
%% 
% Here we'll demonstrate how we can do adaptive beamforming with the USTB.
% As we know from before we can do conventional delay-and-sum.

%% Conventional DAS beamforming is default
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.boxcar;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uff.window.none;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;

b_data_das_CC = bmf.go(postprocess.coherent_compound);
b_data_das_CC.plot(100,['DAS coherent compounded'],80);

b_data_das_IC = bmf.go(postprocess.incoherent_compound);
b_data_das_IC.plot(101,['DAS incoherent compounded'],80);

%%
% If you want to use an adaptive beamforming algorithm you need to make the
% beamforming object an object of one of the subclasses to the *beamformer*
% class. For example the *coherence_factor* subclass.

bmf_cf=coherence_factor();
bmf_cf.channel_data=channel_data;
bmf_cf.scan=sca;

bmf_cf.receive_apodization.window=uff.window.boxcar;
bmf_cf.receive_apodization.f_number=1.7;
bmf_cf.receive_apodization.apex.distance=Inf;

bmf_cf.transmit_apodization.window=uff.window.none;
bmf_cf.transmit_apodization.f_number=1.7;
bmf_cf.transmit_apodization.apex.distance=Inf;

b_data_cf = bmf_cf.go(postprocess.coherent_compound);
b_data_cf.plot(102,['CF'],80);

%%
% However, if you allready have defined yourself a beamformer object you
% can simply call the constructor to the subclass with the precious object
% as input. We'll show this with the *phase_coherence_factor* subclass.
bmf_pcf=phase_coherence_factor(bmf);

% The PCF have a paramter *gamma* that we can set. If we don't set it gets
% a default value.
bmf_pcf.gamma = 1;
b_data_pcf = bmf_pcf.go(postprocess.coherent_compound);
b_data_pcf.plot(103,['PCF'],80);