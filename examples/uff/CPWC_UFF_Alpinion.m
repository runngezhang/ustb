%% Reading CPWC data from an UFF file recorded from an Alpinion scanner
%
% In this example we show how to read channel and beamformed data from a
% UFF (Ultrasound File Format) file recorded with an Alpinion scanner.
% You will need an internet connection to download data.
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net>
%  and Muyinatu Lediju Bell <mledijubell@jhu.edu> 27.05.2017_

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

clear all; close all;
% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data
addpath(local_path);

% We have to different Alpinion CPWC datasets, comment out the one to use
filename='Alpinion_L3-8_CPWC_hyperechoic_scatterers.uff';
%filename='Alpinion_L3-8_CPWC_hypoechoic.uff';

% check if the file is available in the local path & downloads otherwise
tools.download(filename, url, local_path);

%% Reading channel data
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
% If the file had beamformed data, let's read that and the channeldata, 
% before we can jump straight to displaying the data since it's allready
% beamformed

if has_b_data
    b_data=uff_file.read('/b_data');
    channel_data=uff_file.read('/channel_data');
else
    %%
    % If it doesn't have any beamformed data at least it should have some
    % channel_data. So let's read that.
    
    channel_data=uff_file.read('/channel_data');
    
    %%
    %
    % And then do the normal routine of defining the scan,
    
    sca=uff.linear_scan();
    sca.x_axis = linspace(channel_data.probe.x(1),channel_data.probe.x(end),512).';
    sca.z_axis = linspace(1e-3,50e-3,512).';
    
    %%
    %
    % setting up and running the beamforming
    bmf=beamformer();
    bmf.channel_data=channel_data;
    bmf.scan=sca;
    
    bmf.receive_apodization.window=uff.window.tukey25;
    bmf.receive_apodization.f_number=1.7;
    bmf.receive_apodization.apex.distance=Inf;
    
    b_data=bmf.go({process.das_mex process.coherent_compounding});
   
    %%
    %
    % Now we can save this beamformed image to that file, so that we don't
    % have to wait for the beamforming next time.
    uff_file.write(b_data,'b_data');
end

%% Display image
%
% And finally display the image.
b_data.plot([],strrep(filename,'_',' '));

%% Write info about channel data
%
% Let's look at the info given about this dataset
channel_data.print_info()
