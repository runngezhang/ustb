clear all;
close all;

%% Define a reconstruction object
recons=reconstruction();

% define a linear scan 
recons.scan=linear_scan();
recons.scan.x_axis=linspace(-11.1e-3,11.1e-3,256).';               % x vector [m]
recons.scan.z_axis=linspace(20e-3,55e-3,512).';                    % z vector [m]

% define the transmit & receive beams
%F-number, transmit apodization, steering angle [rad], length of the edge smoothing area [elements], order of the edge smoothing polynomial
recons.transmit_beam=beam(1,E.apodization_type.none,0,0);
recons.receive_beam=beam(1.2,E.apodization_type.hanning,0,20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plane wave
load('../../data/joris/lv_neonatal_planewave.mat');

% translating values from FieldSim3 mat file
c0=sim.propagation.c;
x0=(1:sim.probe.N_elements_az).*sim.probe.pitch_az; x0=x0-mean(x0);
geom=[x0.' zeros(sim.probe.N_elements_az,2)];
time=((1:size(sim.data.channel_data,1)).'-1)/sim.fs;
angle=sim.scan.scanSeq{1}.txBeam.tilt(1,1);
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,3) size(sim.data.channel_data,6)]);

% define dataset
cpw_dataset=cpw('Joris data, CPW',...           % dataset name
                 E.signal_format.RF,...         % signal format (RF/IQ)
                 c0,...                         % reference speed of sound (m/s)
                 angle,...                      % vector of plane wave angles (rad)
                 time,...                       % time vector (s)
                 data,...                       % data [samples,channels,firings,frames]
                 geom);                         % probe geometry [x y z] (m)
             
% transform signal format to IQ 
cpw_dataset.demodulate(true);       

% request reconstruction 
cpw_dataset.image_reconstruction(recons);

% write a video to disk
% filename, dynamic_range, video_size_in_pixels 
recons.write_video('joris_iq_cpw_linear_scan.avi',40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% diverging wave
load('../../data/joris/lv_neonatal_divergingwave.mat');

% translating values from FieldSim3 mat file
c0=sim.propagation.c;
x0=(1:sim.probe.N_elements_az).*sim.probe.pitch_az; x0=x0-mean(x0);
geom=[x0.' zeros(sim.probe.N_elements_az,2)];
desired_opening_angle = 60 * pi/180;
virtual_source_depth = (sim.probe.transducerSize_az/2) / tan(desired_opening_angle/2);
source=[0 0 -virtual_source_depth];
time=((1:size(sim.data.channel_data,1)).'-1)/sim.fs+virtual_source_depth/c0; %% < - we must solve the problem of the origin convention
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,3) size(sim.data.channel_data,6)]);

% define dataset
dw_dataset=vs('Joris data, DW',...              % dataset name
                 E.signal_format.RF,...         % signal format (RF/IQ)
                 c0,...                         % reference speed of sound [m/s]
                 source,...                     % vector of virtual source positions [x y z] (m)
                 time,...                       % time vector (s)
                 data,...                       % data [samples,channels,firings,frames]
                 geom);                         % probe geometry [x y z] (m)
             
% transform signal format to IQ 
dw_dataset.demodulate(true); 

% request reconstruction
dw_dataset.image_reconstruction(recons);

% write a video to disk
recons.write_video('joris_iq_dw_linear_scan.avi',40);
