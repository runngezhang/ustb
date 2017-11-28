%% Reading FI data from an UFF file recorded from a Verasonics Scanner
%
% In this example we show how to read channel data from a
% UFF (Ultrasound File Format) file recorded with a Verasonics scanner.
% You will need an internet connection to download data.
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net>
%   and Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%
%   $Last updated: 2017/10/06$

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
filename='L7_FI_Verasonics.uff';

% checks if the data is in your data path, and downloads it otherwise.
% The defaults data path is under USTB's folder, but you can change this
% by setting an environment variable with setenv(DATA_PATH,'the_path_you_want_to_use');
tools.download(filename, url, data_path);   

%% Reading data
%
% Let's first check if we are lucky and the file allready contains
% beamformed_data that we can display.
display=true;
content = uff.index([data_path filesep filename],'/',display);


%% Channel data
% If it doesn't have any beamformed data at least it should have some
% channel_data. So let's read that.

channel_data=uff.read_object([data_path filesep filename],'/channel_data');

%%
figure();plot(channel_data.data(:,1));

%%
% dmod = preprocess.demodulation();
% dmod.input = channel_data;
% channel_data = dmod.go()
%%
%
% And then do the normal routine of defining the scan,
x_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    x_axis(n)=channel_data.sequence(n).source.x;
end
z_axis=linspace(1e-3,55e-3,512).';
scan=uff.linear_scan('x_axis',x_axis,'z_axis',z_axis);

%%
%
% setting up and running the pipeline
mid=midprocess.das();
mid.dimension = dimension.both();

mid.channel_data=channel_data;
mid.scan=scan;

mid.transmit_apodization.window=uff.window.scanline;

mid.receive_apodization.window=uff.window.none;
mid.receive_apodization.f_number=1.7;

b_data=mid.go();

%% Display image
%
% And finally display the image.
b_data.plot([],'Beamformed image');

img = b_data.get_image('none');
figure;imagesc(db(abs(fftshift(fft2(ifftshift(img))))));
title('No MLA no apod');
colormap jet;
%% Beamforming with MLA's
MLA = 2;

scan_MLA=uff.linear_scan('x_axis',linspace(x_axis(1),x_axis(end),length(x_axis)*MLA)','z_axis',z_axis);

mid_MLA=midprocess.das();
mid_MLA.dimension = dimension.receive();

mid_MLA.channel_data=channel_data;
mid_MLA.scan=scan_MLA;

mid_MLA.transmit_apodization.window=uff.window.boxcar;
mid_MLA.transmit_apodization.MLA = MLA;
mid_MLA.transmit_apodization.MLA_overlap = MLA;
mid_MLA.transmit_apodization.f_number = 1.7;
mid_MLA.transmit_apodization.probe = channel_data.probe;
mid_MLA.transmit_apodization.sequence = channel_data.sequence;
mid_MLA.transmit_apodization.focus = scan_MLA;
tx_apo = mid_MLA.transmit_apodization.data;
%%
tx_apo_matrix = reshape(tx_apo,scan_MLA.N_z_axis,scan_MLA.N_x_axis,128);

f999 = figure(999);imagesc(tx_apo_matrix(:,:,55))
title('Tx apodization MLA/RTB');
saveas(f999,'tx_apod_img')
%%
mid_MLA.receive_apodization.window=uff.window.boxcar;
mid_MLA.receive_apodization.f_number=1.7;
mid_MLA.receive_apodization.probe=channel_data.probe;
mid_MLA.receive_apodization.focus=scan_MLA;
rx_apo = mid_MLA.receive_apodization.data;
rx_apo_matrix = reshape(rx_apo,scan_MLA.N_z_axis,scan_MLA.N_x_axis,128);

figure;imagesc(rx_apo_matrix(:,:,50))
%%
b_data_MLA=mid_MLA.go();
b_data_MLA.plot(777,'Beamformed image MLA');


%%
post_sum = postprocess.coherent_compounding();
post_sum.input = b_data_MLA;
b_data_final_image = post_sum.go();
b_data_final_image.plot();


%%
b_data_MLA.write('data/rtb_example_artifact_with_apod.uff','/b_data')
% b_data_MLA.plot([],'Beamformed image MLA',40);
% ylim(gca,[25 35])

%%
img = b_data_MLA.get_image('none');
figure;imagesc(db(abs(fftshift(fft2(ifftshift(img))))));
title('MLA with spline Boxcar');colormap jet
