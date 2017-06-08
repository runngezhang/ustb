%% FI simulation on a linear array with the USTB built-in Fresnel simulator and using Receive Processes
% 
% In this example we show how to use the built-in fresnel simulator in USTB
% to generate a Conventional Focused Imaging (single focal depth) dataset 
% for a linear array and a linear scan and show how it can be beamformed 
% with USTB using several different receive processes that USTB has to offer.
%
% This tutorial assumes familiarity with the contents of the 
% <../../linear_array/html/CPWC_linear_array.html 'CPWC simulation with the 
% USTB built-in Fresnel simulator'> tutorial. Please feel free to refer 
% back to that for more details.
% 
% _by Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no> and Arun
% Asokan Nair <anair8@jhu.edu> 14.03.2017_

%%
% 
% Clear the memory of any lingering settings and data, and close all 
% previously opened plots.

clear all;
close all;

%% Phantom
%
% We start off defining an appropriate *phantom* structure to image. 
% Our phantom here is simply a single point scatterer. USTB's implementation 
% of *phantom* comes with a *plot* method to visualize the phantom for free!

pha=uff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[0,  0, 40e-3, 1];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% Probe
%
% The next UFF structure we look at is *probe*. It contains information 
% about the probe's geometry. USTB's implementation of *probe* comes with a 
% *plot* method too. When combined with the previous figure we can see the
% position of the probe respect to the phantom.

prb=uff.linear_array();
prb.N=128;                  % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=5000e-6; % element height [m]
prb.plot(fig_handle);

%% Pulse
% 
% We then define the pulse-echo signal which is done here using the 
% *fresnel* simulator's *pulse* structure. We could also use 
% <http://field-ii.dk/ 'Field II'> for a more accurate model.

pul=uff.pulse();
pul.center_frequency=5.2e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% Sequence generation
%
% Now, we shall generate our sequence! Keep in mind that the *fresnel* simulator
% takes the same sequence definition as the USTB beamformer. In UFF and
% USTB a sequence is defined as a collection of *wave* structures. 
% 
% For our example here, we define a sequence of 200 focused beams spanning 
% a lateral range of $[-2, 2]$ mm. The focal depth is set as 40 mm. 
% The *wave* structure too has a *plot* method.

N=200;                      % number of focused beams
x_axis=linspace(-2e-3,2e-3,N).';
z0=40e-3;
seq=uff.wave();
for n=1:N 
    seq(n)=uff.wave();
    seq(n).probe=prb;

    seq(n).source.xyz=[x_axis(n) 0 z0];
    
    seq(n).apodization=uff.apodization();
    seq(n).apodization.window=uff.window.tukey50;
    seq(n).apodization.f_number=1.7;
    seq(n).apodization.apex.distance=Inf;
    seq(n).apodization.scan=uff.scan(seq(n).source.xyz);
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% The Fresnel simulator
%
% Finally, we launch the built-in simulator. The simulator takes in a
% *phantom*, *pulse*, *probe* and a sequence of *wave* structures along 
% with the desired sampling frequency, and returns a *channel_data* UFF 
% structure.

sim=fresnel();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();

%% Demodulate the channel data
% 
% USTB has a demodulator class that enables the demodulation of the channel
% data as done below.

dem=demodulator();
dem.channel_data=channel_data;
channel_data=dem.go();
 
%% Scan
%
% The scan area is defines as a collection of pixels spanning our region of 
% interest. For our example here, we use the *sector_scan* structure to 
% generate a sector scan. *scan* too has a useful *plot* method it can call.

z_axis=linspace(39e-3,41e-3,100).';
sca=uff.linear_scan();
for n=1:length(x_axis)
    sca(n)=uff.linear_scan(x_axis(n),z_axis);
    sca(n).plot(fig_handle,'Scenario');    
end
 
%% Beamformer
%
% With *channel_data* and a *scan* we have all we need to produce an
% ultrasound image. We now use a USTB structure *beamformer*, that takes an
% *apodization* structure in addition to the *channel_data* and *scan*.

bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;

%% 
%
% The *beamformer* structure allows you to implement different beamformers 
% by combination of multiple built-in *processes*. By changing the *process*
% chain other beamforming sequences can be implemented. It returns yet 
% another *UFF* structure: *beamformed_data*.
% 
% To achieve the goal of this example, we use delay (implemented in 
% the *delay_matlab()* process) and then stack each line into a single image.

b_data=bmf.go({process.delay_matlab() process.stack()});

%% Test out some of the many receive beamforming processes USTB has to offer

%% 
%
% First, we can use the *coherent_compounding* process to coherently
% compound the data and then display it.
cc=process.coherent_compounding();
cc.beamformed_data=b_data;
cc_data=cc.go();
cc_data.plot([],cc.name);

%% 
%
% Or, we could use the *incoherent_compounding* process to incoherently
% compound the data and then display it.
ic=process.incoherent_compounding();
ic.beamformed_data=b_data;
ic_data=ic.go();
ic_data.plot([],ic.name);

%% 
%
% Or, we could take the *max* process to take the max across the received
% data and then display it.
mv=process.max();
mv.beamformed_data=b_data;
mv_data=mv.go();
mv_data.plot([],mv.name);

%% 
%
% We could also use the *coherence_factor* process which implements the 
% Mallart-Fink coherence factor beamforming to beamform the data.
cf=process.coherence_factor();
cf.channel_data=bmf.channel_data;
cf.transmit_apodization=bmf.transmit_apodization;
cf.receive_apodization=bmf.receive_apodization;
cf.beamformed_data=b_data;
cf_data=cf.go();
cf.CF.plot([],'Mallart-Fink Coherence factor',60,'none'); % show the coherence factor
cf_data.plot([],cf.name);

%% 
%
% Alternatively, we could use the *phase_coherence_factor* process which 
% implements the Camacho-Fritsch phase coherence factor beamforming method.
% We are truly spoilt for choice!
pcf=process.phase_coherence_factor();
pcf.channel_data=bmf.channel_data;
pcf.transmit_apodization=bmf.transmit_apodization;
pcf.receive_apodization=bmf.receive_apodization;
pcf.beamformed_data=b_data;
pcf_data=pcf.go();
pcf.FCC.plot([],'Camacho-Fritsch Phase coherence factor',60,'none'); % show the phase coherence factor
pcf_data.plot([],pcf.name);