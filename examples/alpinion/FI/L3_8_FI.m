% Read FI linear dataset from the Alpinion Research scanner
%
%   This examples shows how to load FI data recorded on the Alpinion
%   Research scanner with L3-8_OMHR_Focused.py example sequence the into a 
%   channel_data object and beamform it with the USTB routines. The data 
%   should be in the path /data/Alpinion/FI_linear/ from the root folder of
%   the USTB. The data can be downloaded from: 
%   https://www.dropbox.com/sh/wv7sa88dvhfszcl/AABYgOJixkKvQM_Q3jJlvFVWa?dl=0
%
%   date : 11.07.2017
%   author: Ole Marius Hoel Rindal (olemarius@olemarius.net)

clear all; close all

local_path = [ustb_path(),'/data/']; % location of example data in this computer

tag = 'FI_hyperechoic_scatterers';  % Dataset 1
%tag = 'FI_hypoechoic';             % Dataset 2
data_folder  = [local_path,'/Alpinion/FI_linear/',tag];

%% initiate Alpinion object pointing to files with data
alp = alpinion();
alp.data_folder = data_folder;

%% Create channel_data object
number_of_frames = 4;
start_frame = 2;
channel_data = alp.read_FI_linear(number_of_frames,start_frame);
%% Define image scan
z_axis=linspace(1e-3,55e-3,512).';
sca=uff.linear_scan();
idx = 1;
for n=1:numel(channel_data.sequence)
    sca(n)=uff.linear_scan('x_axis',channel_data.sequence(n).source.x,'z_axis',z_axis);
end
%% Define BEAMFORMER pipeline
pipe=pipeline();
pipe.channel_data=channel_data;
pipe.scan=sca;
pipe.receive_apodization.window=uff.window.tukey25;
pipe.receive_apodization.f_number=1.7;

% Do beamforming
b_data=bmf.go({midprocess.das_mex postprocess.stack});
%% Show image
b_data.plot(11,['Alpinion FI dataset'],60);
%%
answer = questdlg('Do you want to save this dataset?');
if strcmp(answer,'Yes')
    if strfind(tag, 'scatterers')
        channel_data.name = {'FI dataset of hyperechoic cyst and points scatterers recorded on an Alpinion scanner with a L3-8 Probe from a CIRS General Purpose Ultrasound Phantom'};
    else
        channel_data.name = {'FI dataset of hypoechic cyst recorded on an Alpinion scanner with a L3-8 Probe from a CIRC General Purpose Ultrasound Phantom'};
    end
    channel_data.author = {'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Muyinatu Lediju Bell <mledijubell@jhu.edu>'};
    channel_data.reference = {'www.ultrasoundtoolbox.com'};
    channel_data.version = {'1.0.2'};
    
    %%
    uff_filename = ['./Alpinion_L3-8_',tag,'.uff']
    channel_data.write(uff_filename,'channel_data');
end