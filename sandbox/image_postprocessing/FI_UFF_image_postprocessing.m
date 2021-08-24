clear all;close all;
url = ['https://drive.google.com/uc?export=download' ...
    '&id=1Nf4Zk0lYo0emV-6YdhHJZz9noqtb35Wb'];  % if not found download from here
        %'&id=19OyvPCP4qUiTECFpUe8r3r_Ys2281j0N'];
% Choose dataset
name = 'Verasonics_P2-4_parasternal_long_subject_1.uff';
% Create full filepath
file = fullfile(data_path(), name);
% check if the file is available in the local path or downloads otherwise
tools.download(file, url)

% If issue with the download tool, download manually from, and move to USTB/data: 
% https://drive.google.com/file/d/1Nf4Zk0lYo0emV-6YdhHJZz9noqtb35Wb/view?usp=sharing
% read the data
channel_data = uff.read_object(file,'/channel_data');
%scan = uff.read_object(file,'/scan'); Potentially read the scan from the file

%% Create the images of the heart.
% Only use 5 frames
channel_data.N_frames = 5;

% Create image scan
depth_axis=linspace(0e-3,110e-3,512).';
azimuth_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    azimuth_axis(n) = channel_data.sequence(n).source.azimuth;
end
scan=uff.sector_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);

% Setup beamformer and beamform image
mid=midprocess.das();
mid.dimension = dimension.both;

mid.channel_data=channel_data;
mid.scan=scan;

mid.transmit_apodization.window=uff.window.scanline;
mid.receive_apodization.window=uff.window.tukey25;
mid.receive_apodization.f_number=1.75;
b_data_das = mid.go();

% Display conventional DAS image
b_data_das.plot([],['DAS'],[],[],[],[],[],'dark')

%% Run image postprocessing 
impost = postprocess.image_postprocessing();
impost.input = b_data_das;
impost.median_m = 5;
impost.median_n = 5;
impost.scan = b_data_das.scan;
impost.channel_data = channel_data;
b_data_post_processed = impost.go();
b_data_post_processed.plot([],['img postprocessed'],[],[],[],[],[],'dark')
caxis([-50 10])

%% Alternating view of processed frames
b_data_combined = uff.beamformed_data(b_data_das);
idx = 1;
for f = 1:b_data_das.N_frames
    b_data_combined.data(:,:,:,idx) = b_data_das.data(:,:,:,f);
    b_data_combined.data(:,:,:,idx+1) = b_data_post_processed.data(:,:,:,f);
    idx = idx + 2;
end
b_data_combined.plot([],['Alternating images'],[],[],[],[],[],'dark')
