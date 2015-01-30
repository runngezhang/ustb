clear all;
close all;

%% Create a reconstruction object
recons=reconstruction();

% define the scan -> only linear scan svailable for the moment 
recons.scan.x_axis=linspace(-2e-3,2e-3,200).';                  % x vector [m]
recons.scan.z_axis=linspace(39e-3,41e-3,100).';                 % z vector [m]

% define the transmit & receive beams
%F-number, transmit apodization, steering angle [rad], length of the edge smoothing area [elements], order of the edge smoothing polynomial
recons.transmit_beam=beam(1.75,E.apodization_type.boxcar,0,0,0);
recons.receive_beam=beam(1.75,E.apodization_type.boxcar,0,0,0);

%% Synthetic transmit aperture
% format RF 
load('../../data/ps/ps_sta_rf.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/ps/
s=sta(s.name,s.format,s.c0,s.time,s.data,s.geom);               % define STA dataset object

recons.name='STA, RF, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show(60);                                                % show 

recons.name='STA, RF, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_1=recons.show(60);                                           % show 

% format IQ
load('../../data/ps/ps_sta_iq.mat');                                      % load data; available at http://folk.ntnu.no/alfonsom/data/ps/
s=sta(s.name,s.format,s.c0,s.time,s.data,s.geom,s.modulation_frequency);  % define STA dataset object

recons.name='STA, IQ, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show(60);                                                % show 

recons.name='STA, IQ, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_2=recons.show(60);                                           % show 


%% Coherent plane wave 
% format RF
load('../../data/ps/ps_cpw_rf.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom);       % define CPW dataset object

recons.name='CPW, RF, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show(60);                                                % show 

recons.name='CPW, RF, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_3=recons.show(60);                                           % show 

% format IQ
load('../../data/ps/ps_cpw_iq.mat');                            % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom,s.modulation_frequency); % define CPW dataset object

recons.name='CPW, IQ, Matlab';                                  % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show(60);                                                % show 

recons.name='CPW, IQ, Mex';                                     % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_4=recons.show(60);                                           % show 

%% Virtual Source
% format RF
load('../../data/ps/ps_vs_rf.mat');                             % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom);       % define VS dataset object

recons.name='VS, RF, Matlab';                                   % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show(60);                                                % show 

recons.name='VS, RF, Mex';                                      % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_5=recons.show(60);                                           % show 

% format IQ
load('../../data/ps/ps_vs_iq.mat');                             % load data; available at http://folk.ntnu.no/alfonsom/data/ps
s=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom,s.modulation_frequency); % define VS dataset object

recons.name='VS, IQ, Matlab';                                   % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.simple_matlab);  % request reconstruction
recons.show(60);                                                % show 

recons.name='VS, IQ, Mex';                                      % reconstruction name (optional)
s.image_reconstruction(recons,E.implementation.mex);            % request reconstruction
im_6=recons.show(60);                                           % show 

figure;
subplot(3,2,1); imshow(im_1); caxis([-60 0]);
subplot(3,2,2); imshow(im_2); caxis([-60 0]);
subplot(3,2,3); imshow(im_3); caxis([-60 0]);
subplot(3,2,4); imshow(im_4); caxis([-60 0]);
subplot(3,2,5); imshow(im_5); caxis([-60 0]);
subplot(3,2,6); imshow(im_6); caxis([-60 0]);

