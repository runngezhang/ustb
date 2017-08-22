%% Read channel data from a CPWC dataset from UFF file
channel_data=uff.read_object('./CPWC_PSF.uff','/channel_data');

%% Define scan (the pixels in the image)
scan=uff.linear_scan();
scan.x_axis = linspace(-2e-3,2e-3,200).';
scan.z_axis = linspace(19e-3,21e-3,100).';

%% Set up beamformer for CPWC
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=scan;

bmf.receive_apodization.window=uff.window.boxcar;
bmf.receive_apodization.f_number=1.2;
bmf.receive_apodization.origo.distance=Inf;

bmf.transmit_apodization = bmf.receive_apodization;

% Beamform using delay-and-sum process and coherent compounding postprocess
cpwc_data=bmf.go({process.das_matlab() process.coherent_compounding()});

% Display image
cpwc_data.plot([],'CPWC');