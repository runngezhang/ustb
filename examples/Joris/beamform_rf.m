clear all;
close all;

%% Definition of the reconstruction area
x_vector = linspace(-7.6e-3,7.6e-3,256);                        % x corrdinates of pixels [m]
z_vector = linspace(0,60e-3,1280);                              % z coordinates of pixels [m]
[x z]=meshgrid(x_vector,z_vector);                              % matrix of pixel locations
recons=struct('x',x,'z',z);                                     % reconstruction structure         

%% Definition of the imaging beam
transmit=struct('FN',1,...                                      % transmit F-number 
                'apodization',E.apodization_type.none,...       % transmit apodization
                'steer_angle',0,...                             % transmit steering angle [rad]
                'damping',0,...                                 % length of the damping area [elements]
                'damping_order',0);                             % order of the damping polynomial
            
receive=struct('FN',1.1,...                                     % receive F-number 
                'apodization',E.apodization_type.hanning,...    % receive apodization 
                'steer_angle',0,...                             % transmit steering angle [rad]
                'damping',15,...                                % length of the damping area [elements]
                'damping_order',2);                             % order of the damping polynomial
            
beam=struct('transmit',transmit,'receive',receive);             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plane wave
load('../../data/joris/lv_neonatal_planewave.mat');

% translating values
c0=sim.propagation.c;
x0=(1:sim.probe.N_elements_az).*sim.probe.pitch_az; x0=x0-mean(x0);
geom=[x0.' zeros(sim.probe.N_elements_az,2)];
time=((1:size(sim.data.channel_data,1)).'-1)/sim.fs;
angle=sim.scan.scanSeq{1}.txBeam.tilt(1,1);
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,3) size(sim.data.channel_data,6)]);

% define dataset
cpw_dataset=cpw('Joris data, CPW',...
                 E.signal_format.RF,...
                 c0,...
                 angle,...
                 time,...
                 data,...
                 geom);

% request reconstruction
[sig,im_cpwi]=cpw_dataset.image_reconstruction(beam,recons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% diverging wave
load('../../data/joris/lv_neonatal_divergingwave.mat');

% translating values
c0=sim.propagation.c;
x0=(1:sim.probe.N_elements_az).*sim.probe.pitch_az; x0=x0-mean(x0);
geom=[x0.' zeros(sim.probe.N_elements_az,2)];
desired_opening_angle = 60 * pi/180;
virtual_source_depth = (sim.probe.transducerSize_az/2) / tan(desired_opening_angle/2);
source=[0 0 -virtual_source_depth];
time=((1:size(sim.data.channel_data,1)).'-1)/sim.fs+virtual_source_depth/c0; %% < - we must solve the problem of the origin convention
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,3) size(sim.data.channel_data,6)]);

% define dataset
dw_dataset=vs('Joris data, CPW',...
                 E.signal_format.RF,...
                 c0,...
                 source,...
                 time,...
                 data,...
                 geom);

% request reconstruction
[sig,im_dw]=dw_dataset.image_reconstruction(beam,recons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating videos

% making a video of CPW
writerObj = VideoWriter('cpw.avi');
open(writerObj);

im_dB=20*log10(im_cpwi./max(im_cpwi(:)));
figure; set(gca,'fontsize',16);
for f=1:size(im_cpwi,3)
    pcolor(recons.x*1e3,recons.z*1e3,im_dB(:,:,f)); shading flat; axis equal; axis tight; colormap gray; caxis([-40 0]);colorbar; hold on;
    xlabel('x [mm]');
    ylabel('z [mm]');
    set(gca,'YDir','reverse');
    set(gca,'fontsize',16);
    axis([min(recons.x(:)) max(recons.x(:)) min(recons.z(:)) max(recons.z(:))]*1e3);
    title(sprintf('CPW, Frame %d',f)); 
    set(gcf,'Position',[969 49 944 1068]);
    drawnow;
    
    frame = getframe(gcf,[0 0 944 1068]);
    writeVideo(writerObj,frame);
    
    hold off;
end

close(writerObj);

% making a video of DW
writerObj = VideoWriter('dw.avi');
open(writerObj);

im_dB=20*log10(im_dw./max(im_dw(:)));
figure; set(gca,'fontsize',16);
for f=1:size(im_dw,3)
    pcolor(recons.x*1e3,recons.z*1e3,im_dB(:,:,f)); shading flat; axis equal; axis tight; colormap gray; caxis([-40 0]);colorbar; hold on;
    xlabel('x [mm]');
    ylabel('z [mm]');
    set(gca,'YDir','reverse');
    set(gca,'fontsize',16);
    axis([min(recons.x(:)) max(recons.x(:)) min(recons.z(:)) max(recons.z(:))]*1e3);
    title(sprintf('DW, Frame %d',f)); 
    set(gcf,'Position',[969 49 944 1068]);
    drawnow;
    
    frame = getframe(gcf,[0 0 944 1068]);
    writeVideo(writerObj,frame);
    
    hold off;
end

close(writerObj);

