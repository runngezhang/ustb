%% Acoustic Radiation Force Imaging from UFF file recorded with the Verasonics ARFI_L7 example
%
% In this example we demonstrate the speckle tracking of displacement created
% by a shear wave induced with an acoustic "push" pulse. Also known as
% acoustic radiation force imaging, or share wave elastography. You will 
% need an internet connection to download the data. Otherwise, you can run 
% the *ARFI_L7.m* Verasonics example and create your own .uff file.
%
% _Ole Marius Hoel Rindal <olemarius@olemarius.net> 07.07.2017_

%% Reading the channel data from the UFF file
clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found data will be downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data on this computer
addpath(local_path);
filename='ARFI_dataset.uff';

% check if the file is available in the local path & downloads otherwise
tools.download(filename, url, local_path);

%%
% Reading channel data from the UFF file
% and print out information about the dataset.

uff_file=uff(filename)
channel_data = uff_file.read('/channel_data');
channel_data.print_info

%% Beamform images 
% First, we need to beamform the images from the channel data. We'll do the
% usual drill of defining the scan and the beamformer object.

% SCAN
sca=uff.linear_scan();
sca.x_axis = linspace(channel_data.probe.x(1),channel_data.probe.x(end),256).';
sca.z_axis = linspace(0,30e-3,256).';
 
% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uff.window.tukey50;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;

% Do beamforming
b_data=bmf.go({process.das_mex process.coherent_compounding});

%% Show beamformed images
% The b-mode images is really not that interesting. The displacement
% created by the shear waves are to small to be seen. Notice that we barely
% see a circular structure at about x = 10 and z = 15 mm. This is a sphere
% in the phantom which is somewhat harder than the surrounding structure.

b_data.plot(1,['B-mode'],[30]);

%% Estimate displacement with speckle tracking
% To actually see the shear wave propagation we need to estimate the
% displacement. This is done using a USTB *process* called
% *pulsed_doppler_speckle_tracking*. Have a look at its references to see
% the details.

% Notice that we call the speckle tracking implemenetation
% *pulsed_doppler_speckle_tracking* this is because it's basically the same
% algorithm that is used for color doppler blood flow imaging. Anyway...

pdst = process.pulsed_doppler_speckle_tracking();
pdst.channel_data = channel_data;
pdst.beamformed_data = b_data;
pdst.z_gate = 4;        % Nbr of samples to average estimate in depth / z
pdst.x_gate = 2;        % Nbr of samples to average estimate in lateral / x
pdst.packet_size = 6;   % How many frames to use in the estimate
displacement_estimation = pdst.go();

%% Display the displacement estimate
% Now, we can finally show the estimated displacement. This is nice to
% visualize as a movie.
%
% In the movie we can clearly see the *acoustic radiation force* push
% that creates the share waves. The push is centered at x = 0 mm and 
% z = 14 mm. Notice how the shear waves interacts with the harder
% spherical structure at x = 10 and z = 15 mm. The wavefront is moving
% faster through the harder tissue and some complex wavefronts are created.
%
% Notice that we give the figure and the title as arguments to the plot
% function. The empty brackets [] is because we don't want to specify any
% dynamic range (well do that with the caxis function). And the 'none' is
% because we don't want to do any compression of the displacement data

f2 = figure(2);clf;
handle = displacement_estimation.plot(f2,['Displacement'],[],'none');
caxis([-0.1*10^-6 0.2*10^-6]); % Updating the colorbar
colormap(gca(f2),'hot');       % Changing the colormap

%% We can do it all in once!!
% Since the *pulsed_doppler_speckle_tracking* is a process, we can trust
% the default paramters (which are the same as the ones used above) and do
% the beamforming and the displacement estimation all in one call!
%
% Isn't the USTB great?!

disp = bmf.go({process.das_mex process.coherent_compounding ...
                                process.pulsed_doppler_speckle_tracking});

%%
% Display the displacement 
% Which gives us the same result as above.
f3 = figure(3);clf;
disp.plot(f3,['Displacement'],[],'none');
caxis([0.01*10^-6 0.3*10^-6]); % Updating the colorbar
colormap(gca(f3),'hot');       % Changing the colormap
