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
N=10;                           % number of diverging waves
x0=linspace(-10e-3,10e-3,N);
for n=1:N 
    seq(n)=wave();
    seq(n).probe=prb;
    seq(n).source.xyz=[x0(n) 0 -10e-3];
    
    %seq(n).apodizator.apodization_window=window.hamming;
    %seq(n).apodizator.origin=seq(n).source;
    %seq(n).apodizator.scan=scan(x0(n),0,20e-3);
    %seq(n).apodizator.f_number=1;
    
    seq(n).sound_speed=pha.sound_speed;
    %seq(n).source.plot(fig_handle,'Scenario');
    %seq(n).plot(); 
    %pause();
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
 
% %check how does it look
% for n=1:N 
%     raw.plot(n);
%     pause();
% end

%% DEMODULATOR
dem=demodulator();
dem.raw_data=raw;
dem_raw=dem.go();
dem_raw=raw;

%% SCAN
%
% This is the generic scan class. Handling will be simplified with children classes for
% linear_scan, sector_scan, volumetric_scan, or so on.
sca=scan();
%x_axis=linspace(min(prb.x),max(prb.x),256);
x_axis=linspace(-5e-3,5e-3,256);
z_axis=linspace(15e-3,25e-3,256);
[X Z]=meshgrid(x_axis,z_axis);

sca.x=X(:);
sca.y=0.*X(:);
sca.z=Z(:);

sca.plot(fig_handle,'Scenario');    % show mesh
 
%% BEAMFORMER
%
% First approximation to the general beamformer
w0=2*pi*dem_raw.modulation_frequency;

%% beamforming
sig=zeros(sca.N_pixels,numel(seq));
wb=waitbar(0,'Beamforming');
for ntx=1:numel(seq)
    waitbar(ntx/numel(seq));
    TF=sqrt((seq(ntx).source.x-sca.x).^2+(seq(ntx).source.y-sca.y).^2+(seq(ntx).source.z-sca.z).^2);
    if ~isinf(seq(ntx).source.distance)
        if(seq(ntx).source.z<0)
            TF=TF-seq(ntx).source.distance;
        else
            TF=TF+seq(ntx).source.distance;
        end
    end
    for nrx=1:prb.N_elements
        RF=sqrt((prb.x(nrx)-sca.x).^2+(prb.y(nrx)-sca.y).^2+(prb.z(nrx)-sca.z).^2);
        delay=(RF+TF)/dem_raw.sound_speed;
        phase_shift=-exp(1i.*w0*delay);
        sig(:,ntx)=sig(:,ntx)+phase_shift.*interp1(dem_raw.time,dem_raw.data(:,nrx,ntx),delay,'linear',0);
    end
    
%     TF0=sqrt((seq(ntx).source.x-0).^2+(seq(ntx).source.y-0).^2+(seq(ntx).source.z-20e-3).^2);
%     if ~isinf(seq(ntx).source.distance)
%         if(seq(ntx).source.z<0)
%             TF0=TF0-seq(ntx).source.distance;
%         else
%             TF0=seq(ntx).source.distance+TF0;
%         end
%     end
%     RF0=sqrt((prb.x-0).^2+(prb.y-0).^2+(prb.z-20e-3).^2);
%     
%     figure(111);
%     subplot(1,2,1)
%     imagesc(1:prb.N_elements,raw.time,raw.data(:,:,ntx)); hold on; grid on;
%     plot(1:prb.N_elements,(RF0+TF0)/raw.sound_speed,'r-');
%     title(ntx);
%     axis tight;
%     subplot(1,2,2)
%     imagesc(x_axis*1e3,z_axis*1e3,abs(reshape(sig(:,ntx),[length(z_axis) length(x_axis)]))); axis equal tight;
%     
%     pause();
    
end
close(wb);

beamformed_data=sum(reshape(sig,[length(z_axis) length(x_axis) size(sig,2)]),3);

%beamformed_data=reshape(sig(:,1),[length(z_axis) length(x_axis) ]);

% convert to intensity values
if(dem_raw.modulation_frequency>0)
    envelope=abs(beamformed_data);
else
    envelope=abs(hilbert(beamformed_data));
end
envelope_dB=20*log10(envelope./max(envelope(:)));

figure3 = figure('Color',[1 1 1]); 
imagesc(x_axis*1e3,z_axis*1e3,envelope_dB); axis tight equal; 
box('on'); 
xlabel('x [mm]');
ylabel('z [mm]')
set(figure3,'InvertHardcopy','off');
caxis([-60 0]); colorbar; colormap gray;
set(gca,'color','black')
%ylim([0 40]);

