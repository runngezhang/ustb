%% Reading data from an UFF file recorded with the Verasonics CPWC_L7 example
%
% In this example we show how to read channel and beamformed data from a
% UFF (Ultrasound File Format) file recorded with the Verasonics example.
% You will need an internet connectionto download data. Otherwise, you can
% run the *CPWC_L7.m* Verasonics example so the file 'L7_CPWC_193328.uff'
% is in the current path.
%
% _by Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no> 
%  and Ole Marius Hoel Rindal <olemarius@olemarius.net>_ 
%
% Last update: 07.08.2017

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

% data location
url='http://ustb.no/datasets/';      % if not found data will be downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data in this computer
filename='L7_CPWC_193328.uff';

% check if the file is available in the local path & downloads otherwise
tools.download(filename, url, local_path);

%% Checking what's inside
%
% Now that the file is in the machine we can start loading data. The first 
% would be to check what is in there with the *uff.index* function 
uff.index([local_path filename],'/',display);

% 
% We see there is a *beamformed_data* dataset with name _b_data_. Let us
% load it and display it. 

b_data=uff.beamformed_data();
b_data.read([local_path filename],'/b_data');
b_data.plot();

%%
%
%
% Let's read the channel data now,
    
channel_data=uff.read_object([local_path filename],'/channel_data');

%%
%
% define a scan
   
sca=uff.linear_scan();
sca.x_axis = linspace(channel_data.probe.x(1),channel_data.probe.x(end),256).';
sca.z_axis = linspace(0,50e-3,256).';
    
%%
%
% and beamform
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
    
bmf.receive_apodization.window=uff.window.none;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.origo=uff.point('xyz',[0 0 -Inf]);
    
bmf.transmit_apodization.window=uff.window.none;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.origo=uff.point('xyz',[0 0 -Inf]);
    
b_data2=bmf.go({process.delay_mex process.coherent_compounding});
b_data2.plot();
