%% Load a CPWI dataset from Fieldsim data and beamforms it with USTB
%
% This example shows how to load the data from a FieldSim simulation into a
% CPWI class, demodulate and beamform it with the USTB routines.  The
% FielsSim toolbox must be installed in the system. To run this example the 
% following data file is needed 
%
%           http://users.ugent.be/~jehvcauw/lv_neonatal_planewave.mat
%
% which can be downloaded from the author's webpage. 
%
% date:     11.03.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%           Joris Van Cauwenberge <Joris.VanCauwenberge@ugent.be>

load('../../data/joris/lv_neonatal_planewave.mat');

%% converting FielsSim data to USTB format

% time vector
time=(0:(size(sim.data.channel_data,1)-1)).'/sim.fs;

% angle vector
angles=zeros(length(sim.scan.scanSeq),1);
for kk = 1:length(sim.scan.scanSeq)
    angles(kk)=sim.scan.scanSeq{kk}.txBeam.tilt(1,1);
end

% data
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,4) size(sim.data.channel_data,6)]);

%% define the CPW dataset and demodulate
cpw_dataset=cpw('Neonatal, CPW',...                       % dataset name
                 E.signal_format.RF,...                   % signal format (RF/IQ)
                 sim.propagation.c,...                    % reference speed of sound (m/s)
                 angles,...                               % vector of plane wave angles (rad)
                 time,...                                 % time vector (s)
                 data,...                                 % data [samples,channels,firings,frames]
                 squeeze(sim.probe.getElementCenters())); % probe geometry [x y z] (m)

% transform signal format to IQ 
cpw_dataset.demodulate(true);   

%% Define a reconstruction object
recons=reconstruction();

% define a linear scan 
recons.scan=linear_scan();
recons.scan.x_axis=linspace(-11.1e-3,11.1e-3,256).';               % x vector [m]
recons.scan.z_axis=linspace(20e-3,55e-3,512).';                    % z vector [m]

% define the synthetic orientation
F_number=1.2;
recons.orientation=orientation();
recons.orientation.transmit_beam=beam(1,E.apodization_type.none);
recons.orientation.receive_beam=beam(F_number,E.apodization_type.hamming,0,20);

% request reconstruction 
cpw_dataset.image_reconstruction(recons);

% write a video to disk
% filename, dynamic_range, video_size_in_pixels 
recons.write_video('neonatal_cpw_linear_scan.avi',40);
