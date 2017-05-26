%% Reading data from an UFF file recorded from a Alpinion scanner
%
% In this example we show how to read channel and beamformed data from a
% UFF (Ultrasound File Format) file recorded with an Alpinion scanner.
% You will need an internet connection to download data.
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net>
%  and Muyinatu Lediju Bell <mledijubell@jhu.edu> 26.05.2017_

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

% data location
url='http://ustb.no/datasets/';      % if not found data will be downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data in this computer
addpath(local_path);

% We have to different Alpinion CPWC datasets, comment out the one to use
filename='Alpinion_L3-8_FI_hyperechoic_scatterers.uff';
%filename='Alpinion_L3-8_FI_hypoechoic.uff';

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
%
%
% If it doesn't have any beamformed data at least it should have some
% channel_data. So let's read that.
% If the file did have beamformed data, will will jump straight to the
% display part later on :)
if ~has_b_data
    
    % Read the channel data
    channel_data=uff_file.read('/channel_data');
    
    %%
    %
    % And then the normal routine of defining the scan,
    
    z_axis=linspace(1e-3,55e-3,512).';
    sca=uff.linear_scan();
    idx = 1;
    for n=1:numel(channel_data.sequence)
        sca(n)=uff.linear_scan(channel_data.sequence(n).source.x,z_axis);
    end
    
    %%
    %
    % setting up and running the beamforming
    bmf=beamformer();
    bmf.channel_data=channel_data;
    bmf.scan=sca;
    
    bmf.receive_apodization.window=uff.window.tukey25;
    bmf.receive_apodization.f_number=1.7;
    bmf.receive_apodization.apex.distance=Inf;
    
    b_data=bmf.go({process.das_mex process.stack});
    %%
    %
    % Now we can save this beamformed image to that file, so that we don't have
    % to wait for the beamforming again.
    uff_file.write(b_data,'b_data');
end

%%
%
% And finally display the image.
b_data.plot();

%%
%
% Let's look at the info given about this dataset
channel_data.print_info()
