%% Example of beamforming with USTB
%
% Example on how to launch beamform data with USTB. Only mex implementations 
% are launched. To run this example the following data files are needed
%
%       ps_sta_rf.mat
%       ps_sta_iq.mat
%       ps_cpw_rf.mat
%       ps_cpw_iq.mat
%       ps_vs_rf.mat
%       ps_vs_iq.mat
%
% which must be located within MATLAB's path.
%
% date:     11.03.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%           Andreas �stvik <andreas.ostvik@ntnu.no>

clear all; close all; clear classes; clc;

%% Create a reconstruction object
recons=reconstruction();

% define the scan 
recons.scan=linear_scan();
recons.scan.x_axis=linspace(-2e-3,2e-3,200).';                  % x vector [m]
recons.scan.z_axis=linspace(39e-3,41e-3,100).';                 % z vector [m]

% define the transmit & receive beams
F_number=1.75;
recons.orientation=orientation();
recons.orientation.transmit_beam=beam(F_number, pyus.apodization_type.boxcar);
recons.orientation.receive_beam=beam(F_number, pyus.apodization_type.boxcar);


%% Synthetic transmit aperture
% format IQ
load('ps_sta_iq.mat');                                      % load data; available at http://folk.ntnu.no/alfonsom/data/ps/
s=sta(s.name,s.format,s.c0,s.time,s.data,s.geom,s.modulation_frequency);  % define STA dataset object

recons.name='STA, IQ, Matlab';                                  % reconstruction name (optional)
pyus.image_reconstruction('ps_sta_iq.mat', recons, pyus.implementation.pyus, 'sta'); % request reconstruction
im_1=recons.show();                                                % show 

%% Coherent plane wave 
% format IQ
load('ps_cpw_iq.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/ps/
s=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom,s.modulation_frequency); % define CPW dataset object
recons.name='CPW, IQ, PyUS'; % reconstruction name (optional)   

pyus.image_reconstruction('ps_cpw_iq.mat', recons, pyus.implementation.pyus, 'cpw');
im_2=recons.show();  % show 

%% Virtual Source
% format IQ
load('ps_vs_iq.mat');                             % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom,s.modulation_frequency); % define VS dataset object

recons.name='VS, IQ, PyUS';                                   % reconstruction name (optional)
pyus.image_reconstruction('ps_vs_iq.mat', recons, pyus.implementation.pyus, 'vs');  % request reconstruction
im_3 = recons.show();                                                % show 

