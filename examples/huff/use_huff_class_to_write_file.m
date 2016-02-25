%% Using HUFF class to write 
%
% This example shows how to interact with the huff class to store data in a 
% hdf5 file according to the specification HUFF v0.0.2. To run this example 
% the following data file is needed
%
%       ps_sta_iq.mat
%
% which must be located within MATLAB's path.

% date:     11.03.2015
% authors:  Alfonso Rodriguez-Molares (alfonsom@ntnu.no)

%% create a huff object
recording=huff('recording.h5','w');

%% Create a reconstruction object
sta_recons=reconstruction('STA, IQ, Mex');

% define the scan -> only linear scan svailable for the moment 
sta_recons.scan=linear_scan();
sta_recons.scan.x_axis=linspace(-2e-3,2e-3,200).';                  % x vector [m]
sta_recons.scan.z_axis=linspace(39e-3,41e-3,100).';                 % z vector [m]

% define the transmit & receive beams
sta_recons.orientation=orientation();
sta_recons.orientation.transmit_beam=beam(1.75,E.apodization_type.boxcar);
sta_recons.orientation.receive_beam=beam(1.75,E.apodization_type.boxcar);

%% Synthetic transmit aperture
% format IQ
load('../../data/ps/ps_sta_iq.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/data/ps/
sta_dataset=sta(s.name,s.format,s.c0,s.time,s.data,s.geom,s.modulation_frequency);  % define STA dataset object

% reconstruction
sta_dataset.image_reconstruction(sta_recons);                           
sta_recons.show();                                                  

% adding datasets to staging area
recording.stage(sta_dataset); 
recording.stage(sta_recons);

%% Coherent Plane Wave
% format IQ
load('../../data/ps/ps_cpw_iq.mat');                              % load data; available at http://folk.ntnu.no/alfonsom/data/ps/
cpw_dataset=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom,s.modulation_frequency);  % define CPW dataset object

cpw_recons = reconstruction('CPW, IQ, Mex',sta_recons);           % new reconstruction, copy template from previous
cpw_dataset.image_reconstruction(cpw_recons);                     % request reconstruction
cpw_recons.show();                                                % show 

% adding datasets to staging area
recording.stage(cpw_dataset);
recording.stage(cpw_recons);

%% Virtual source
% format RF
load('../../data/ps/ps_vs_rf.mat');                                 % load data; available at http://folk.ntnu.no/alfonsom/data/ps/
vs_dataset=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom);  % define VS dataset object

vs_recons=reconstruction('VS, RF, Mex',sta_recons);                 % new reconstruction, copy template from previous
vs_dataset.image_reconstruction(vs_recons);                         % request reconstruction
vs_recons.show();                                                   % show 

% adding datasets to staging area
recording.stage(vs_dataset);
recording.stage(vs_recons);

%% writting hdf5 to disk
recording.write();


