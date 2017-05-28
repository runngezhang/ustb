%% PICMUS challenge: simulation, resolution-distortion test
%
% This example reads (or downloads if it cannot find the data locally) a 
% dataset used in the <http://ieeexplore.ieee.org/document/7728908/ PICMUS challenge>
% and beamforms it with USTB's general beamformer.
% A 75 plane-wave sequence was simulated with <http://field-ii.dk/ Field
% II> to estimate the beamforming method resolution and geometric
% distortion. 
%
% _by Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no> 
%  and Olivier Bernard <olivier.bernard@insa-lyon.fr> 26.05.2017_

%% Checking if the file is in the local path, and downloading otherwise
%
% We define the local path and the url where the data is stored

% data location
url='http://ustb.no/datasets/';      % if not found data will be downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data in this computer
addpath(local_path);
filename='PICMUS_simulation_resolution_distortion.uff';

% check if the file is available in the local path & downloads otherwise
tools.download(filename, url, local_path);

%% Reading channel data
%
% Now that the file is on the path we deine a *UFF* object to interact
% with it.

uff_file=uff(filename,'read');

%%
%
% This dataset should contain the following structures:
% * *channel_data*,
% * *beamformed_data* and,
% * *scan*
%
% We can check it out with the *index* function
display=true;
content = uff_file.index('/',display);

%%
%
% This dataset should contain the following structures:
% * *channel_data*,
% * *beamformed_data* and,
% * *scan*
%
% We can check it out with the *index* function
display=true;
content = uff_file.index('/',display);

%% Plotting beamformed_data
%
% We can read the *beamformed_data* object and plot it 

b_data=uff_file.read('/beamformed_data');
b_data.plot();

%% Loading channel data & scan
%
% The file also contain channel_data and scan. We read it so we can
% replicate the beamformed image in the UFF file.

channel_data=uff_file.read('/channel_data');
scan=uff_file.read('/scan');

%% Beamforming
%
% We define a beamformer, and the corresponding transmit and apodization
% windows, and launch it.

bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=scan;
    
% receive apodization
bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;

% transmit apodization
bmf.transmit_apodization.window=uff.window.tukey50;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;

% launch beamforming
b_data_new=bmf.go({process.das_mex process.coherent_compounding});

%% Comparing results
%
% We plot both images side by side.

figure;
b_data.plot(subplot(1,2,1),'Original');
b_data_new.plot(subplot(1,2,2),'New');

