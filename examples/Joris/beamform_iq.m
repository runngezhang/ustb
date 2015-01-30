clear all;
close all;

%% Define a reconstruction object
recons=reconstruction();

% define the scan -> only linear scan svailable for the moment 
recons.scan.x_axis=linspace(-7.6e-3,7.6e-3,256).';               % x vector [m]
recons.scan.z_axis=linspace(0,60e-3,512).';                      % z vector [m]

% define the transmit & receive beams
%F-number, transmit apodization, steering angle [rad], length of the edge smoothing area [elements], order of the edge smoothing polynomial
recons.transmit_beam=beam(1,E.apodization_type.none,0,0,0);
recons.receive_beam=beam(1.1,E.apodization_type.hanning,0,15,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plane wave
load('../../data/joris/lv_neonatal_planewave.mat');

% translating values
c0=sim.propagation.c;
x0=(1:sim.probe.N_elements_az).*sim.probe.pitch_az; x0=x0-mean(x0);
geom=[x0.' zeros(sim.probe.N_elements_az,2)];
time=((1:size(sim.data.channel_data,1)).'-1)/sim.fs;
angle=sim.scan.scanSeq{1}.txBeam.tilt(1,1);
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,3) size(sim.data.channel_data,6)]);

% define dataset
cpw_dataset=cpw('Joris data, CPW',...
                 E.signal_format.RF,...
                 c0,...
                 angle,...
                 time,...
                 data,...
                 geom);
             
% transform signal format to IQ 
cpw_dataset.demodulate(true);       

% request reconstruction 
cpw_dataset.image_reconstruction(recons);

% write a video to disk
% filename, dynamic_range, video_size_in_pixels 
recons.write_video('joris_iq_cpw.avi',40,[500 1000]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% diverging wave
load('../../data/joris/lv_neonatal_divergingwave.mat');

% translating values
c0=sim.propagation.c;
x0=(1:sim.probe.N_elements_az).*sim.probe.pitch_az; x0=x0-mean(x0);
geom=[x0.' zeros(sim.probe.N_elements_az,2)];
desired_opening_angle = 60 * pi/180;
virtual_source_depth = (sim.probe.transducerSize_az/2) / tan(desired_opening_angle/2);
source=[0 0 -virtual_source_depth];
time=((1:size(sim.data.channel_data,1)).'-1)/sim.fs+virtual_source_depth/c0; %% < - we must solve the problem of the origin convention
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,3) size(sim.data.channel_data,6)]);

% define dataset
dw_dataset=vs('Joris data, VS',...
                 E.signal_format.RF,...
                 c0,...
                 source,...
                 time,...
                 data,...
                 geom);
             
% transform signal format to IQ 
dw_dataset.demodulate(true); 

% request reconstruction
dw_dataset.image_reconstruction(recons);

% write a video to disk
recons.write_video('joris_iq_dw.avi',40,[500 1000]);
