%% Reading data from an UFF file recorded with the Verasonics CPWC_L7 example
%
% In this example we show how to read channel and beamformed data from a
% UFF (Ultrasound File Format) file recorded with the Verasonics example.
% You will need an internet connectionto download data. Otherwise, you can
% run the *CPWC_L7.m* Verasonics example so the file 'L7_CPWC_193328.uff'
% is in the current path.
%
% _by Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no> 15.05.2017
%  and Ole Marius Hoel Rindal <olemarius@olemarius.net>_

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

% data location
url='http://ustb.no/datasets/';   % if not found data will be downloaded from here
local_path='./';                  % location of example data in this computer
filename='L7_CPWC_193328.uff';

% check if the file is available in the local path & downloads otherwise
tools.download(filename, url, local_path);

%% Reading channel data
%
% Now that the file is on the path let us create a *UFF* object to interact
% with the file. We open it in "append" mode.


uff_file=uff(filename);


%%
%
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
    
    channel_data=uff_file.read('/channel_data');
    
    %%
    %
    % And then the normal routine of defining the scan,
    
    sca=uff.linear_scan();
    sca.x_axis = linspace(channel_data.probe.x(1),channel_data.probe.x(end),256).';
    sca.z_axis = linspace(0,50e-3,256).';
    
    %%
    %
    % setting up and running the beamforming
    bmf=beamformer();
    bmf.channel_data=channel_data;
    bmf.scan=sca;
    
    bmf.receive_apodization.window=uff.window.none;
    bmf.receive_apodization.f_number=1.7;
    bmf.receive_apodization.apex.distance=Inf;
    
    bmf.transmit_apodization.window=uff.window.none;
    bmf.transmit_apodization.f_number=1.7;
    bmf.transmit_apodization.apex.distance=Inf;
    
    b_data=bmf.go({process.delay_matlab process.coherent_compounding});
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


