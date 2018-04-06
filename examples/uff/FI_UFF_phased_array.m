clear all;close all;
url='http://ustb.no/datasets/';      % if not found downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data
addpath(local_path);

% Choose dataset
filename='Verasonics_P2-4_parasternal_long.uff';
% check if the file is available in the local path or downloads otherwise
tools.download(filename, url, local_path);
channel_data = uff.read_object([local_path, filename],'/channel_data');
scan = uff.read_object([local_path, filename],'/scan');

%%
% Print info about the dataset. Remeber that if you want to use this dataset
% you have to reference this article!
channel_data.print_authorship

%% Do beamforming

channel_data.N_frames = 1;

depth_axis=linspace(0e-3,110e-3,512).';
azimuth_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    azimuth_axis(n) = channel_data.sequence(n).source.azimuth;
end

scan=uff.sector_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);

mid=midprocess.das();
mid.dimension = dimension.both();
mid.channel_data=channel_data;
mid.scan=scan;

mid.transmit_apodization.window=uff.window.scanline;
mid.receive_apodization.window=uff.window.tukey25;
mid.receive_apodization.f_number = 1.7;
% This will result in a beamformed_data object with the delayed and not
% summed channel data.
b_data = mid.go();

%% Plot the two images
b_data.plot(3,['DAS']);
mid.receive_apodization.plot([],32);