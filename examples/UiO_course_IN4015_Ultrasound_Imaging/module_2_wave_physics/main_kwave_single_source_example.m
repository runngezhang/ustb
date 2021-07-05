% Receive beamforming on k-wave single source exampler
%
% This example demonstrates how to run a simulation in k-wave recording a
% single originating from a single source. We are then showing how we can
% import the recorded signal into a channel data UFF object in the USTB and
% beamforming the data into an image. Thus we will do "receive beamforming"
% and reconstruct an image of the single source. 
%
% Your task is to implement your own pixelbased receive beamformer. For
% background you can read page 1 and 2 of Jørgen Grythes document
% "Beamforming Algorithms - beamformers" or pages 22-29. However, remember
% that here you only need to to receive beamforming.
%
% Author: Ole Marius Hoel Rindal
clearvars;

% =========================================================================
% Run K-wave SIMULATION
% =========================================================================
% Set up some parameters
number_of_sensors = 4; % Define how many receive sensors with lambda/2 spacing
dynamic_range = 40; % How many decibels to display in image
%transmit_signal = 'sinus';
transmit_signal = 'gaussian_pulse';

[channel_data, kgrid] = run_kwave_simulation(number_of_sensors,transmit_signal)

%%  
% Reconstructing an image of the single source using the USTB.

% Create scan UFF object to define the set of pixels we want to beamform
scan = uff.linear_scan();
scan.x_axis = linspace(kgrid.x_vec(1),kgrid.x_vec(end),1024)';
scan.z_axis = linspace(0,30e-3,1024)';

% Perform the beamforming using the DAS midprocess object
das = midprocess.das();
das.channel_data = channel_data;
das.scan = scan;
b_data = das.go()

% Visualise the image
b_data.plot([],[],[dynamic_range])

%% Your own receive beamforming
% Now, your assignment is to implement a receive beamformer. However, most
% of the code is allready written, so you simply have to get the receive
% delay correct (thus finish line 67) and your image should be similar to
% the one resulting from the USTB.

% Allow me to extract the variables you need
x_element_position = channel_data.probe.x; %The x-position of the elements
z_element_position = channel_data.probe.z; %The z-position of the elements

x_pixels = reshape(scan.x, scan.N_x_axis, scan.N_z_axis); %The x position of the pixels
z_pixels = reshape(scan.z', scan.N_x_axis, scan.N_z_axis); %The z position of the pixels

ch_data = hilbert(channel_data.data); %The raw channel data in analytical form

% Empty variables of correct dimension of variables you need to calculate
receive_delay = zeros(scan.N_x_axis,scan.N_z_axis,channel_data.N_elements);
delayed_data = zeros(scan.N_x_axis,scan.N_z_axis,channel_data.N_elements);
img = zeros(scan.N_x_axis,scan.N_z_axis);

for rx = 1:channel_data.N_elements
    % See equation (2) in Grythe, or equation (1.15) in Rindal (remember to convert (1.15) to seconds).
    receive_delay(:,:,rx) = 0;
    delayed_data(:,:,rx) = interp1(channel_data.time,ch_data(:,rx),receive_delay(:,:,rx),'linear',0);
    img = img + delayed_data(:,:,rx);
end

%% Plotting the image from the USTB and the resulting images from your beamformer
figure;
b_data.plot(subplot(121),['USTB image'],[dynamic_range])
subplot(122)
imagesc(scan.x_axis*1000,scan.z_axis*1000,db(abs(img./max(img(:)))))
xlabel('x [mm]');ylabel('z [mm]')
colormap gray; caxis([-dynamic_range 0])
colorbar; axis image;
title('Your image');
set(gca,'fontsize',14);

%% More indepth analysis of the partial and final results
% To be commented...
figure()
subplot(211)
imagesc(scan.x_axis*1000,scan.z_axis*1000,x_pixels);title('x-pixels');
xlabel('x [mm]');ylabel('z [mm]');axis image;
subplot(212)
imagesc(scan.x_axis*1000,scan.z_axis*1000,z_pixels);title('z-pixels');
xlabel('x [mm]');ylabel('z [mm]');axis image;

% Plot the delayed signal from each sensor
figure(10);
for i = 1:number_of_sensors
    subplot(number_of_sensors/2,number_of_sensors/2,i)
    imagesc(scan.x_axis*1000,scan.z_axis*1000,real(delayed_data(:,:,i)));
    xlabel('x [mm]');ylabel('z [mm]');
    title(['Delayed signal from sensor ',num2str(i)])
end

figure;
imagesc(scan.x_axis*1000,scan.z_axis*1000,sum(real(delayed_data),3))
xlabel('x [mm]');ylabel('z [mm]')
figure;
imagesc(scan.x_axis*1000,scan.z_axis*1000,db(abs(img./max(img(:)))))
xlabel('x [mm]');ylabel('z [mm]')
colormap gray; caxis([-20 0])