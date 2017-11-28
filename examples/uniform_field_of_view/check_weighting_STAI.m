

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
filename='FieldII_speckle_simulation.uff';

% checks if the data is in your data path, and downloads it otherwise.
% The defaults data path is under USTB's folder, but you can change this 
% by setting an environment variable with setenv(DATA_PATH,'the_path_you_want_to_use');
tools.download(filename, url, data_path);   


%% Channel data
% If it doesn't have any beamformed data at least it should have some
% channel_data. So let's read that.

channel_data=uff.read_object([data_path filesep filename],'/channel_data');
scan=uff.linear_scan('x_axis',linspace(channel_data.probe.x(1),channel_data.probe.x(end),256).','z_axis',linspace(1e-3,15e-3,256).');


f_numbers = [2 5];
windows = [uff.window.boxcar uff.window.tukey25 uff.window.hamming]
for f_number = f_numbers
    for window = windows
        
        %Setting all channel data to ones
        channel_data.data = ones(size(channel_data.data));
        
        % setting up and running the pipeline
        mid=midprocess.das();
        mid.dimension = dimension.both();
        
        mid.channel_data=channel_data;
        mid.scan=scan;
        
        mid.transmit_apodization.window=window;
        mid.transmit_apodization.f_number=1.7;
        
        mid.receive_apodization.window=window;
        mid.receive_apodization.f_number=1.7;
        
        b_data_das = mid.go();
        

        %% Calculate weighting using a function under tools
        [apod] = tools.uniform_fov_weighting(mid);
        das_image = b_data_das.get_image('none-complex');

        
        figure
        subplot(211)
        imagesc(real(das_image));colorbar
        title(sprintf('Image %s, no weighting, F number=%d',mid.receive_apodization.window,mid.receive_apodization.f_number(1)));
        subplot(212)
        imagesc(abs(das_image./apod));colorbar
        title(sprintf('Image %s, weighting, F number=%d',mid.receive_apodization.window,mid.receive_apodization.f_number(1)));
        %saveas(f101,['Figures/FindWeighting/',sprintf('img_%s_F_number=%d',mid.receive_apodization.window,mid.receive_apodization.f_number(1))],'png');

    end
end