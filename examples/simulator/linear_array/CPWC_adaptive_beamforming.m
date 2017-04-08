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

N=5;                           % number of plane waves
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
 
%% Beamformer
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=0;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uff.window.tukey50;
bmf.transmit_apodization.f_number=0;
bmf.transmit_apodization.apex.distance=Inf;

%% 
% Here we'll demonstrate how we can do adaptive beamforming with the USTB.
% As we know from before we can do conventional delay-and-sum by using the
% default matlab implementation as a parameter to the *go* method of the
% *beamformer* structure as shown below.

%% Conventional DAS
figure(1);
axis_handle = subplot(131);

b_data_das=bmf.go(@bmf.matlab,@postprocess.coherent_compound);

% show
b_data_das.plot(axis_handle,['Delay-and-sum (DAS)'],80);


%%
% If you want to use an adaptive beamforming algorithm you give a oject to
% an adaptive beamforming implementation as a paramtere to the *go* method.
% The adaptive beamforming implementation have to be a subclass of the
% *ADAPTIVE_BEAMFORMER* class. These are all found under the namespace
% folder *adaptive_beamformers*. Below we show how to use both the
% *coherence_factor* and *phase_coherence_factor* implementation in USTB.
%% CF beamforming
axis_handle = subplot(132);

b_data_CF=bmf.go(adaptive_beamformers.coherence_factor(),@postprocess.coherent_compound);

% show
b_data_CF.plot(axis_handle,['Coherence Factor (CF)'],80);

%% PCF beamforming
% The phase coherence factor has an parameter *gamma* that the user can
% set. 
axis_handle = subplot(133);
pcf = adaptive_beamformers.phase_coherence_factor();
pcf.gamma = 1;
b_data_pcf=bmf.go(pcf,@postprocess.coherent_compound);

% show
b_data_pcf.plot(axis_handle,['Phase Coherence Factor (PCF)'],80);