%% Reading CPWC data from an UFF file recorded from an Alpinion scanner
%
% In this example we show how to read channel and beamformed data from a
% UFF (Ultrasound File Format) file recorded with an Alpinion scanner.
% You will need an internet connection to download data.
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net>
%  and Muyinatu Lediju Bell <mledijubell@jhu.edu>_
%
% Last update 07.08.2017

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data

% We have to different Alpinion CPWC datasets, comment out the one to use
short_filename='Alpinion_L3-8_CPWC_hyperechoic_scatterers.uff';
filename=[local_path short_filename];

% check if the file is available in the local path & downloads otherwise
tools.download(short_filename, url, local_path);

%% Reading channel data
%
% Let's first check if we are lucky and the file allready contains
% beamformed_data that we can display.
display=true;
content = uff.index(filename,'/',display);

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
    b_data=uff.read_object(filename,'/b_data');
    b_data.plot([],'Stored image');
end

%% Beamform data
% If it doesn't have any beamformed data at least it should have some
% channel_data. So let's read that.

channel_data=uff.read_object(filename,'/channel_data');

%%
%
% And then do the normal routine of defining the scan,

sca=uff.linear_scan();
sca.x_axis = linspace(channel_data.probe.x(1),channel_data.probe.x(end),512).';
sca.z_axis = linspace(1e-3,50e-3,512).';

%%
%
% setting up pipeline
pipe=pipeline();
pipe.channel_data=channel_data;
pipe.scan=sca;
pipe.receive_apodization.window=uff.window.tukey25;
pipe.receive_apodization.f_number=1.7;

b_data2=pipe.go({midprocess.das_mex postprocess.coherent_compounding});

%% Display image
%
% And finally display the image.
b_data2.plot([],'Beamformed image');

%% Write info about channel data
%
% Let's look at the info given about this dataset
channel_data.print_authorship();
