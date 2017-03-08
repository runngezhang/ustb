% Script to develop simulator & general beamformer

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2017/02/23$

clear all;
%close all;

%% PHANTOM
pha=phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[0,  0, 20e-3, 1];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% PROBE
%
% This is the generic probe class. Handling will be simplified with children classes for
% linear_array, curvilinear_array, 2D_matrix, etc.
prb=probe();
Nx=128;                                 % number of elements in each row
Ny=1;                                   % number of rows
az_pitch=300e-6;                        % probe pitch in azimuth [m]
el_pitch=1500e-6;                       % probe pitch in elevation [m]
az_w=270e-6;                            % element width [m]
%el_w=1470e-6;                           % element height [m]
el_w=5000e-6;                           % element height [m]

[X,Y] = meshgrid((1:Nx)*az_pitch,(1:Ny)*el_pitch);

x0=X(:);                    % element position in the x_axis (m)
y0=Y(:);                    % element position in the y_axis [m]
z0=zeros(Nx*Ny,1);          % element position in the z_axis [m]

x0=x0-mean(x0);
y0=y0-mean(y0);
z0=z0-mean(z0);

theta=zeros(Nx*Ny,1);       % element orientation in the azimuth direction [rad]
phi=zeros(Nx*Ny,1);         % element orientation in the elevation direction [rad]

prb.geometry=[x0 y0 z0 theta phi az_w*ones(Nx*Ny,1) el_w*ones(Nx*Ny,1)]; % probe geometry
prb.plot(fig_handle);

%% PULSE
pul=pulse();
pul.center_frequency=5.2e6;                           % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;                         % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SEQUENCE GENERATION
N=10;                           % number of focussed waves
x0=linspace(-10e-3,10e-3,N);
for n=1:N 
    seq(n)=wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[x0(n) 0 10e-3];
    
    seq(n).apodization.apodization_window=window.hamming;
    seq(n).apodization.apex.distance=Inf;
    seq(n).apodization.f_number=1;
    seq(n).apodization.scan=scan(seq(n).source.x,seq(n).source.y,seq(n).source.z);
    
    seq(n).apodization.apodization_window
    
    seq(n).sound_speed=pha.sound_speed;
%     seq(n).source.plot(fig_handle,'Scenario');
%     seq(n).plot(); 
%     pause();
end

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% we launch the simulation
raw=sim.go();
 
%% DEMODULATOR
dem=demodulator();
dem.raw_data=raw;
dem_raw=dem.go();

%% SCAN
%
sca=linear_scan();
sca.x_axis=linspace(-4e-3,4e-3,256).';
sca.z_axis=linspace(16e-3,24e-3,256).';

sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
%
% First approximation to the general beamformer

bmf=beamformer();
bmf.raw_data=dem_raw;
bmf.scan=sca;
bmf.probe=prb;
bmf.sequence=seq;
bmf.receive_apodization.apodization_window=window.hanning;
bmf.receive_apodization.f_number=1;
bmf.receive_apodization.apex.distance=Inf;

% beamforming
bmf_data=bmf.go();

% show
bmf_data(6).plot();

