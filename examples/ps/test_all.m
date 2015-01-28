clear all;
close all;

%% Definition of reconstruction domain
x_vector=linspace(-2e-3,2e-3,200);      % x vector [m]
z_vector=linspace(39e-3,41e-3,100);     % z vector [m]
[x z]=meshgrid(x_vector,z_vector);      % matrix of the reconstruction locations
recons=struct('x',x,'z',z);             % reconstruction structure         

%% Definition of the imaging beam
transmit=struct('FN',1.75,...                                   % transmit F-number 
                'apodization',E.apodization_type.boxcar,...     % transmit apodization
                'steer_angle',0,...                             % transmit steering angle [rad]
                'damping',15,...                                % length of the damping area [elements]
                'damping_order',2);                             % order of the damping polynomial
            
receive=struct('FN',1.75,...                                    % receive F-number 
                'apodization',E.apodization_type.boxcar,...     % receive apodization 
                'steer_angle',0,...                             % transmit steering angle [rad]
                'damping',15,...                                % length of the damping area [elements]
                'damping_order',2);                             % order of the damping polynomialbeam=struct('transmit',transmit,'receive',receive);             % receive steering angle [rad]

beam=struct('transmit',transmit,'receive',receive);                         

%% Synthetic transmit aperture
% format RF 
load('../../data/ps/ps_sta_rf.mat'); % http://folk.ntnu.no/alfonsom/data/ps/ps_sta_rf.mat
s=sta(s.name,s.format,s.c0,s.time,s.data,s.geom);
[sig,im]=s.image_reconstruction(beam,recons); s.show(recons,im,60);
[sig,im]=s.image_reconstruction(beam,recons,E.implementation.simple_matlab); s.show(recons,im,60);

% format IQ
load('../../data/ps/ps_sta_iq.mat'); % http://folk.ntnu.no/alfonsom/data/ps/ps_sta_iq.mat
s=sta(s.name,s.format,s.c0,s.time,s.data,s.geom,s.modulation_frequency);
[sig,im]=s.image_reconstruction(beam,recons); s.show(recons,im,60);
[sig,im]=s.image_reconstruction(beam,recons,E.implementation.simple_matlab); s.show(recons,im,60);

%% Coherent plane wave 
% format RF
load('../../data/ps/ps_cpw_rf.mat'); % http://folk.ntnu.no/alfonsom/data/ps/ps_cpw_rf.mat
s=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom);
[sig,im]=s.image_reconstruction(beam,recons); s.show(recons,im,60);
[sig,im]=s.image_reconstruction(beam,recons,E.implementation.simple_matlab); s.show(recons,im,60);

% format IQ
load('../../data/ps/ps_cpw_iq.mat'); % http://folk.ntnu.no/alfonsom/data/ps/ps_cpw_iq.mat
s=cpw(s.name,s.format,s.c0,s.angle,s.time,s.data,s.geom,s.modulation_frequency);
[sig,im]=s.image_reconstruction(beam,recons); s.show(recons,im,60);
[sig,im]=s.image_reconstruction(beam,recons,E.implementation.simple_matlab); s.show(recons,im,60);

%% Virtual Source
% format RF
load('../../data/ps/ps_vs_rf.mat'); % http://folk.ntnu.no/alfonsom/data/ps/ps_vs_rf.mat
s=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom);
[sig,im]=s.image_reconstruction(beam,recons); s.show(recons,im,60);
[sig,im]=s.image_reconstruction(beam,recons,E.implementation.simple_matlab); s.show(recons,im,60);

% format IQ
load('../../data/ps/ps_vs_iq.mat'); % http://folk.ntnu.no/alfonsom/data/ps/ps_vs_iq.mat
s=vs(s.name,s.format,s.c0,s.source,s.time,s.data,s.geom,s.modulation_frequency);
[sig,im]=s.image_reconstruction(beam,recons); s.show(recons,im,60);
[sig,im]=s.image_reconstruction(beam,recons,E.implementation.simple_matlab); s.show(recons,im,60);

