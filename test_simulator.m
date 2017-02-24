% Script to develop simulator & general beamformer

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2017/02/23$

clear all;
close all;

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
N=10;                           % number of waves
x=linspace(-10e-3,10e-3,N);
for n=1:N 
    seq(n)=wave();
    seq(n).probe=prb;
    seq(n).source=source();
    seq(n).source.xyz=[x(n) 0 -10e-3];
    seq(n).sound_speed=pha.sound_speed;
    seq(n).source.plot(fig_handle,'Scenario');
end
seq(1).plot(); % plot one of the delay profiles

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;                   % phantom
sim.pulse=pul;                % transmitted pulse
sim.probe=prb;                    % probe
sim.sequence=seq;             % beam sequence
sim.sampling_frequency=41.6e6;     % sampling frequency [Hz]

% we launch the simulation
raw=sim.go();
 
% check how does it look
% for n=1:N 
%     raw.plot(n);
%     pause();
% end

%% SCAN
%
% This is the generic scan class. Handling will be simplified with children classes for
% linear_scan, sector_scan, volumetric_scan, or so on.
sca=scan();
x_axis=linspace(min(prb.x),max(prb.x),128);
z_axis=linspace(0e-3,40e-3,128);
[X Z]=meshgrid(x_axis,z_axis);

sca.x=X(:);
sca.y=0.*X(:);
sca.z=Z(:);

sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
%
% First approximation to the general beamformer
w0=0;

%% beamforming
sig=zeros(sca.N_pixels,numel(seq));
wb=waitbar(0,'Beamforming');
for ntx=1:numel(seq)
    waitbar(ntx/numel(seq));
    TF=sqrt((seq(ntx).source.x-sca.x).^2+(seq(ntx).source.y-sca.y).^2+(seq(ntx).source.z-sca.z).^2)-seq(ntx).source.distance;
    for nrx=1:prb.N_elements
        RF=sqrt((prb.x(nrx)-sca.x).^2+(prb.y(nrx)-sca.y).^2+(prb.z(nrx)-sca.z).^2);
        delay=(RF+TF)/raw.sound_speed;
        phase_shift=exp(1i.*w0*delay);
        sig(:,ntx)=phase_shift.*interp1(raw.time,raw.data(:,nrx,ntx),delay,'linear',0);
    end
end
close(wb);



% convert to intensity values
envelope_drf=abs(hilbert(sr_image));
envelope_drf_dB=20*log10(envelope_drf./max(envelope_drf(:)));

figure3 = figure('Color',[1 1 1]); 
imagesc(x_axis*1e3,raw.time*pha.sound_speed/2*1e3,envelope_drf_dB); axis tight equal; 
box('on'); 
xlabel('x [mm]');
ylabel('z [mm]')
set(figure3,'InvertHardcopy','off');
caxis([-60 0]); colorbar; colormap gray;
set(gca,'color','black')
ylim([0 40]);

