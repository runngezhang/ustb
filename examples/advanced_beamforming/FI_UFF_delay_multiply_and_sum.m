%% Delay Multiply And Sum on FI data from an UFF file
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net> 28.05.2017_

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

clear all; close all;
tic
% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data
addpath(local_path);

% We have two different Alpinion CPWC datasets, comment out the one to use
filename='Alpinion_L3-8_FI_hyperechoic_scatterers.uff';
%filename='Alpinion_L3-8_FI_hypoechoic.uff';

% check if the file is available in the local path or downloads otherwise
tools.download(filename, url, local_path);

%% Reading data
%
% Now that the file is on the path let us create a *UFF* object to interact
% with the file. We open it in "append" mode.

uff_file=uff(filename)

%%
%
% Let's first check if we are lucky and the file allready contains
% beamformed_data that we can display.
display=true;
content = uff_file.index('/',display);

has_b_data = false;
for i = 1:length(content)
    if strcmp(content{i}.class,'uff.beamformed_data')
        has_b_data = true; % We found a beamformed data object!
    end
end

%%
% If it doesn't have any beamformed data at least it should have some
% channel_data. So let's read that.

channel_data=uff_file.read('/channel_data');

%Actually, to save some time, let's hack the channel data and only use one
%of the frames in the dataset.

channel_data.N_frames = 1;

%% Define Scan
%
% And then do the normal routine of defining the scan,

z_axis=linspace(25e-3,45e-3,1024).';
sca=uff.linear_scan();
idx = 1;
for n=1:numel(channel_data.sequence)
    sca(n)=uff.linear_scan(channel_data.sequence(n).source.x,z_axis);
end

%% Set up delay part of beamforming
%
% setting up and running the beamforming
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.boxcar;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uff.window.none;

b_data=bmf.go({process.delay_mex process.stack});

%%
dmas = process.delay_multiply_and_sum()
dmas.dimension = dimension.receive();
dmas.receive_apodization = bmf.receive_apodization;
dmas.transmit_apodization = bmf.transmit_apodization;
dmas.beamformed_data = b_data;
dmas.channel_data = channel_data;
dmas.scan = sca;
b_data_dmas = dmas.go();    % Launch beamformer
b_data_dmas.plot(100,'DMAS'); % Display image

%% DAS
das = process.coherent_compounding();
das.beamformed_data = b_data;
b_data_das = das.go();
b_data_das.plot(2,'DAS');

%% 
% Plot both in same plot with connected axes, try to zoom!
figure(3);
b_data_dmas.plot(subplot(2,1,1),'DMAS'); % Display image
ax(1) = gca;
b_data_das.plot(subplot(2,1,2),'DAS'); % Display image
ax(2) = gca;
linkaxes(ax);

%% Compare resolution
% Plot the lateral line through some of the scatterers

% Let's get the images as a N_z_axis x N_x_axis image
dmas_img = b_data_dmas.get_image();
das_img = b_data_das.get_image();

%%
% So that we can plot the line through the group of scatterers
line_idx = 250;
figure(4);clf;
plot(b_data_dmas.scan.x_axis*10^3,dmas_img(line_idx,:),...
                               'DisplayName','DMAS','LineWidth',2);hold on;
plot(b_data_das.scan.x_axis*10^3,das_img(line_idx,:),...
                               'DisplayName','DAS','LineWidth',2);
xlabel('x [mm]');xlim([0 20]);ylabel('Amplitude [dB]');legend show
title(sprintf('Lateral line through %.2f mm',...
                                  b_data_dmas.scan.z_axis(line_idx)*10^3));
%%
%So that we can plot the line through the signle scatterer
line_idx = 747;
figure(5);clf;
plot(b_data_dmas.scan.x_axis*10^3,dmas_img(line_idx,:),...
                                'DisplayName','DMAS','LineWidth',2);hold on;
plot(b_data_das.scan.x_axis*10^3,das_img(line_idx,:),...
                                'DisplayName','DAS','LineWidth',2);
xlabel('x [mm]');xlim([-15 5]);ylabel('Amplitude [dB]');legend show
title(sprintf('Lateral line through %.2f mm',...
                                  b_data_dmas.scan.z_axis(line_idx)*10^3));
toc();