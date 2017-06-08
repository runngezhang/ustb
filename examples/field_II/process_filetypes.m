clear all; close all;

filetype = 'uff'

if strcmp(filetype,'mat')
    load 'FieldII_PSF_simulation.mat' 'channel_data'
elseif strcmp(filetype,'uff')
    %% Save UFF dataset
    uff_file=uff('FieldII_PSF_simulation.uff');
    channel_data = uff_file.read('/channel_data');
end

%% SCAN
sca=uff.linear_scan(linspace(-4e-3,4e-3,256).', linspace(16e-3,24e-3,256).');

%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;
bmf.receive_apodization.window=uff.window.boxcar;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uff.window.boxcar;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;

% Delay and sum on receive, then coherent compounding
b_data=bmf.go({process.das_mex() process.coherent_compounding()});
% Display imagexs
b_data.plot([],['From ',filetype,' file'])

