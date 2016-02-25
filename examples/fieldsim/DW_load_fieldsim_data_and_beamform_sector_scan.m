%% Load a DW dataset from Fieldsim data and beamforms it with USTB
%
% This example shows how to load the data from a FieldSim simulation into a
% VSI class, demodulate and beamform it with the USTB routines.  The
% FielsSim toolbox must be installed in the system. To run this example the 
% following data file is needed 
%
%           http://users.ugent.be/~jehvcauw/lv_neonatal_divergingwave.mat
%
% which can be downloaded from the author's webpage. 
%
% date:     11.03.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%           Joris Van Cauwenberge <Joris.VanCauwenberge@ugent.be>
clear all;
close all;

load('../../data/joris/lv_neonatal_divergingwave.mat');

%% converting FielsSim data to USTB format

% source location
desired_opening_angle = 60 * pi/180;
virtual_source_depth = (sim.probe.transducerSize_az/2) / tan(desired_opening_angle/2);
source=[0 0 -virtual_source_depth];

% time vector
time=((1:size(sim.data.channel_data,1)).'-1)/sim.fs+virtual_source_depth/sim.propagation.c; %% < - we must solve the problem of the origin convention

% data
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,4) size(sim.data.channel_data,6)]);

%% define the CPW dataset and demodulate
dw_dataset=vs('Neonatal, DW',...                       % dataset name
                 E.signal_format.RF,...                   % signal format (RF/IQ)
                 sim.propagation.c,...                    % reference speed of sound (m/s)
                 source,...                     % vector of virtual source positions [x y z] (m)
                 time,...                                 % time vector (s)
                 data,...                                 % data [samples,channels,firings,frames]
                 squeeze(sim.probe.getElementCenters())); % probe geometry [x y z] (m)

% transform signal format to IQ 
dw_dataset.demodulate(true);   

%% Define a reconstruction object
recons=reconstruction();

% define a linear scan 
recons.scan=sector_scan();
recons.scan.apex=[0 0 0];
recons.scan.azimuth_axis=linspace(-0.3,0.3,256).';      % azimuth vector [rad]
recons.scan.depth_axis=linspace(20e-3,55e-3,512).';     % depth vector [m]    

% define the transmit & receive beams
F_number=1.1;
recons.orientation=orientation();
recons.orientation.transmit_beam=beam(0,E.apodization_type.none,0,0);
recons.orientation.receive_beam=beam(F_number,E.apodization_type.hanning,0,20);

% request reconstruction 
dw_dataset.image_reconstruction(recons);

% write a video to disk
% filename, dynamic_range, video_size_in_pixels 
recons.write_video('neonatal_dw_sector_scan.avi',40,[1280 800],10);
