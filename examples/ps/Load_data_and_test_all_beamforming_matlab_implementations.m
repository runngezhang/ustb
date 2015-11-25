%% Example of beamforming with USTB
%
% Example on how to launch beamform data with USTB. To run this 
% example the following data files are needed
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

%% Create a reconstruction object
recons=reconstruction();

% define the scan 
recons.scan=linear_scan();
recons.scan.x_axis=linspace(-2e-3,2e-3,200).';                  % x vector [m]
recons.scan.z_axis=linspace(39e-3,41e-3,100).';                 % z vector [m]

% define the transmit & receive beams
F_number=1.75;
recons.orientation=orientation();
recons.orientation.transmit_beam=beam(F_number,E.apodization_type.boxcar);
recons.orientation.receive_beam=beam(F_number,E.apodization_type.boxcar);

%% Synthetic transmit aperture
% format RF 
load('../../data/ps/ps_sta_rf.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/ps/
s=sta(s.name,s.format,s.c0,s.time,s.data,s.geom);               % define STA dataset object

recons.name='STA, RF, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show();                                                % show 

recons.name='STA, RF, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_1=recons.show();                                           % show 

% format IQ
load('../../data/ps/ps_sta_iq.mat');                                      % load data; available at http://folk.ntnu.no/alfonsom/data/ps/
s=sta(s.name,s.format,s.c0,s.time,s.data,s.geom,s.modulation_frequency);  % define STA dataset object

recons.name='STA, IQ, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show();                                                % show 

recons.name='STA, IQ, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_2=recons.show();                                           % show 


%% Coherent plane wave 
% format RF
load('../../data/ps/ps_cpw_rf.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom);       % define CPW dataset object

recons.name='CPW, RF, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show();                                                % show 

recons.name='CPW, RF, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_3=recons.show();                                           % show 

% format IQ
load('G:\begrenset\bmt\ustb_data\v1.9\psf\ps_cpw_iq.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom,s.modulation_frequency); % define CPW dataset object

recons.name='CPW, IQ, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show();                                                % show 

recons.name='CPW, IQ, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_4=recons.show();                                           % show 

%% Virtual Source
% format RF
load('../../data/ps/ps_vs_rf.mat');                             % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom);       % define VS dataset object

recons.name='VS, RF, Matlab';                                   % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show();                                                % show 

recons.name='VS, RF, Mex';                                      % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_5=recons.show();                                           % show 

% format IQ
load('../../data/ps/ps_vs_iq.mat');                             % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom,s.modulation_frequency); % define VS dataset object

recons.name='VS, IQ, Matlab';                                   % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show();                                                % show 

recons.name='VS, IQ, Mex';                                      % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_6=recons.show();                                           % show 

figure;
subplot(3,2,1); imshow(im_1); caxis([-60 0]);
subplot(3,2,2); imshow(im_2); caxis([-60 0]);
subplot(3,2,3); imshow(im_3); caxis([-60 0]);
subplot(3,2,4); imshow(im_4); caxis([-60 0]);
subplot(3,2,5); imshow(im_5); caxis([-60 0]);
subplot(3,2,6); imshow(im_6); caxis([-60 0]);

