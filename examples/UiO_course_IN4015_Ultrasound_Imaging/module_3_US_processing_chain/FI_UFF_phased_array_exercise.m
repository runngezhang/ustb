%% Phased array beamforming
% In this example we see how one can do phased array imaging with the USTB.
% In the second part of the script is an exercise where one can implement
% phased array beambased beamforming to get the same image as was produced
% with the USTB.
%
% Author: Ole Marius Hoel Rindal <olemarius@olemarius.net>
% Date: 28.05.21

% Clear up
close all;

% Read the data, poentitally download it
url='http://ustb.no/datasets/';      % if not found downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data
addpath(local_path);

% Choose dataset
filename='Verasonics_P2-4_parasternal_long_small.uff';
%filename='FieldII_P4_point_scatterers.uff';
% check if the file is available in the local path or downloads otherwise
tools.download(filename, url, local_path);
channel_data = uff.read_object([local_path, filename],'/channel_data');

%% Do phased array beamforming with the USTB
% Define the scan
depth_axis=linspace(0e-3,110e-3,1024).';
azimuth_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    azimuth_axis(n) = channel_data.sequence(n).source.azimuth;
end
scan=uff.sector_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);

% Call the midprocessor delay and sum 
mid=midprocess.das();
mid.channel_data=channel_data;
mid.scan=scan;
mid.dimension = dimension.both()
mid.transmit_apodization.window=uff.window.scanline;
mid.receive_apodization.window=uff.window.none;
b_data = mid.go();

%% Plot image
b_data.plot(3,['DAS'],[],[],[],[],[],'dark');

% Get the first frame from UFF beamformed data object
ustb_img = b_data.get_image('none');
ustb_img = ustb_img(:,:,1)./max(max((ustb_img(:,:,1))));

figure; 
imagesc(rad2deg(scan.azimuth_axis), scan.depth_axis*1000, db(abs(ustb_img)))
ylabel('Depth [mm]');xlabel(['angle'])
colormap gray; caxis([-60 0]); title('USTB image in "beamspace"');


%% Implement beambased beamforming from scratch
% Please see the slides on beambased beamforming, and especially the slide
% on the geometry og beambased beamforming on tips on how to implement
% this.
% 
% A hint is to review the exercise from module 2 on wave physics where you
% implemented an receive beamformer. Now you are extending it to also
% include to compensate for the transmit part of the propagation delay as
% well as handling multiple transmits.

% First let us extract the variables you need
N_depth = scan.N_depth_axis;                     % Number of depth pixels
N_transmits = channel_data.N_waves;              % Number of transmits 
N_elements = channel_data.N_elements;            % Number of elements
x_element_position = channel_data.probe.x;       % The x-position of the elements
z_element_position = channel_data.probe.z;       % The z-position of the elements
c = channel_data.sound_speed;                    % Speed of sound
fs = channel_data.sampling_frequency;            % The sampling frequency
rfData = hilbert(channel_data.data(:,:,:,1));    % The RF data as the analytical signal
sample_time = channel_data.time;                 % The time in seconds for each sample

for t = 1:N_transmits 
    angles(t) = channel_data.sequence(t).source.azimuth; % Transmit angles (in radians)
    offset(t) = channel_data.sequence(t).delay;         % Extract the offset for each transmit event
                                        % this offset is used to get
                                        % the correct time zero convention
end

% Let us define the radial distance
radial_distance = linspace(0,110e-3,N_depth);    % Depth (or radial distance)

% Let us define some empty variables to contain the values you calculate
delays = zeros(N_transmits, N_depth, N_elements);% Buffer for the calculated delays
img = zeros(N_depth,N_transmits);                % Buffer for the final image

% Pseudocode for beambased phased array beamformer
% for each transmit
    % calculate the transmit part and the receive part of the delay for
    % every receive channel. Remember to calculate the delays in seconds
    % not distance, and remember to subtract the offset for each transmit
    % event to get a correct time zero convention.

    %for each receive channel
        % beamform by interpolating (using interp1) into the calculated
        % delays for each radial dept sample for each transmit
        % the call to interp1 might look like this:
        % interp1(sample_time, rfData(:,r,t)', delays(t,:,r))'
            % where r is the current receive channel and t is the current transmit
            % you allready have sample_time and rfData, but need to calculate
            % the propriate delays

%Normalize the image to maximum value = 1
img = img./max(img(:));

% Plot the image you created and compare with the USTB image. The
% difference between them should be 0
figure;
subplot(131)
imagesc(db(abs(img)));
ax(1) = gca();
caxis([-60 0])
title('Your image')

subplot(132)
imagesc(db(abs(ustb_img)));
caxis([-60 0])
ax(2) = gca();
title('USTB image')

subplot(133)
imagesc(db(abs(img))-db(abs(ustb_img)));
colorbar
ax(3) = gca();
title('The difference');
linkaxes(ax);

% Finally, let's use the scan convert tool in the USTB to scan convert the
% image from beam space to pixel space.
[img_sc, Xs, Zs] = tools.scan_convert(db(abs(img./max(img(:)))),angles,radius,512,512);

figure;
imagesc(Xs*1000, Zs*1000, img_sc);
ylabel(['z [mm]']); xlabel('x [mm]');
colormap gray;caxis([-60 0])
title('Your image scan converted')
