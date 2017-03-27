% Generating the images for IUS 2017 abstract with USTB built-in Fresnel simulator

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/16$

clear all;
close all;

z0=20e-3;           % depth of the scatterer [m]
Fs=41.6e6;          % sampling frequency [Hz]
F_number=1.2;       % F-number
uwindow=huff.window.boxcar;

%% PHANTOM
pha=huff.phantom();
pha.sound_speed=1540;         % speed of sound [m/s]
pha.points=[0,  0, z0, 1];    % point scatterer position [m]
fig_handle=pha.plot();             
             
%% PROBE
prb=huff.linear_array();
prb.N=128;                  % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=5000e-6; % element height [m]
prb.plot(fig_handle);

%% PULSE
pul=huff.pulse();
pul.center_frequency=5.2e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% SCAN
sca=huff.linear_scan(linspace(-2e-3,2e-3,200).', linspace(z0-1e-3,z0+1e-3,100).');
sca.plot(fig_handle,'Scenario');    % show mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEQUENCE GENERATION
N=128;                      % number of waves
seq=huff.wave();
for n=1:N 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[prb.x(n) prb.y(n) prb.z(n)];
    
    seq(n).apodization.window=huff.window.sta;
    seq(n).apodization.apex=seq(n).source;
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;            % phantom
sim.pulse=pul;              % transmitted pulse
sim.probe=prb;              % probe
sim.sequence=seq;           % beam sequence
sim.sampling_frequency=Fs;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
  
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
bmf.receive_apodization.window=uwindow;
bmf.receive_apodization.f_number=F_number;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uwindow;
bmf.transmit_apodization.f_number=F_number;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
stai_data=bmf.go(@bmf.matlab,@postprocess.coherent_compound);

% show
stai_data.plot([],'STAI');

data=reshape(stai_data.data,[sca.N_z_axis sca.N_x_axis]);
lat_prof=abs(data(50,100:200)); lat_prof=lat_prof./max(lat_prof);
lat_prof=20*log10(lat_prof);
xlat=interp1(lat_prof,sca.x_axis(100:200),-6);
fwhm=2*xlat

figure(101);
plot(sca.x_axis(100:200)*1e3,lat_prof,'linewidth',2); grid on; hold on; axis tight;
%plot([-xlat xlat]*1e3,[-6 -6],'ro')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEQUENCE GENERATION
seq=huff.wave();
for n=1:sca.N_x_axis 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[sca.x_axis(n) 0 z0];
    
    seq(n).apodization.window=uwindow;
    seq(n).apodization.f_number=F_number;
    seq(n).apodization.apex.distance=Inf;
    seq(n).apodization.scan.xyz=seq(n).source.xyz;
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=Fs;      % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% SCAN
fi_sca=huff.linear_scan();
for n=1:sca.N_x_axis
    fi_sca(n)=huff.linear_scan(sca.x_axis(n),sca.z_axis);
    fi_sca(n).plot(fig_handle,'Scenario');    
end
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=fi_sca;
bmf.receive_apodization.window=uwindow;
bmf.receive_apodization.f_number=F_number;
bmf.receive_apodization.apex.distance=Inf;

% beamforming
fi_data=bmf.go(@bmf.matlab,@postprocess.stack);

% show
fi_data.plot([],'FI');

data=reshape(fi_data.data,[sca.N_z_axis sca.N_x_axis]);
lat_prof=abs(data(50,100:200)); lat_prof=lat_prof./max(lat_prof);
lat_prof=20*log10(lat_prof);
xlat=interp1(lat_prof,sca.x_axis(100:200),-6);
fwhm=2*xlat

figure(101);
plot(sca.x_axis(100:200)*1e3,lat_prof,'--','linewidth',2); grid on; hold on;
%plot([-xlat xlat]*1e3,[-6 -6],'ro')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CPWC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEQUENCE GENERATION
angle_max=atan(1/2/F_number);
lambda=pha.sound_speed/pul.center_frequency;
N=round(prb.pitch*prb.N_elements/lambda/F_number);
angles=linspace(-angle_max,angle_max,N);
seq=huff.wave();
for n=1:N 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.azimuth=angles(n);
    seq(n).source.distance=Inf;
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=Fs;      % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uwindow;
bmf.receive_apodization.f_number=F_number;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uwindow;
bmf.transmit_apodization.f_number=F_number;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
cpwc_data=bmf.go(@bmf.matlab,@postprocess.coherent_compound);

% show
cpwc_data.plot([],'CPWC');

data=reshape(cpwc_data.data,[sca.N_z_axis sca.N_x_axis]);
lat_prof=abs(data(50,100:200)); lat_prof=lat_prof./max(lat_prof);
lat_prof=20*log10(lat_prof);
xlat=interp1(lat_prof,sca.x_axis(100:200),-6);
fwhm=2*xlat

figure(101);
plot(sca.x_axis(100:200)*1e3,lat_prof,':','linewidth',2); grid on; hold on;
%plot([xlat]*1e3,[-6],'ro')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DWI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEQUENCE GENERATION
zs=-10e-3;
xs=(z0-zs)*tan(angles);
seq=huff.wave();
for n=1:N 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[xs(n) 0 zs];
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=Fs;      % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uwindow;
bmf.receive_apodization.f_number=F_number;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uwindow;
bmf.transmit_apodization.f_number=F_number;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
dwi_data=bmf.go(@bmf.matlab,@postprocess.coherent_compound);

% show
dwi_data.plot([],'DWI');

data=reshape(dwi_data.data,[sca.N_z_axis sca.N_x_axis]);
lat_prof=abs(data(50,100:200)); lat_prof=lat_prof./max(lat_prof);
lat_prof=20*log10(lat_prof);
xlat=interp1(lat_prof,sca.x_axis(100:200),-6);
fwhm=2*xlat

figure(101);
plot(sca.x_axis(100:200)*1e3,lat_prof,'.-','linewidth',2); grid on; hold on;
%plot([xlat]*1e3,[-6],'ro')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RTB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEQUENCE GENERATION
zs=40e-3;
xs=(zs-z0)*tan(angles);
seq=huff.wave();
for n=1:N 
    seq(n)=huff.wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[xs(n) 0 zs];
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% SIMULATOR
sim=simulator();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=Fs;      % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uwindow;
bmf.receive_apodization.f_number=F_number;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uwindow;
bmf.transmit_apodization.f_number=F_number;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
rtb_data=bmf.go(@bmf.matlab,@postprocess.coherent_compound);

% show
rtb_data.plot([],'RTB');

data=reshape(rtb_data.data,[sca.N_z_axis sca.N_x_axis]);
lat_prof=abs(data(50,100:200)); lat_prof=lat_prof./max(lat_prof);
lat_prof=20*log10(lat_prof);
xlat=interp1(lat_prof,sca.x_axis(100:200),-6);
fwhm=2*xlat

figure(101);
plot(sca.x_axis(100:200)*1e3,lat_prof,'--','linewidth',2); grid on; hold on;
legend('STAI','FI','CPWC','DWI','RTB');
xlabel('x [mm]');
ylabel('Intensity [dB]');
set(gca,'fontsize', 14);
title('Lateral profile');
