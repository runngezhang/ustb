%% Read CPWI dataset from the Alpinion Research scanner
%
%   This examples shows how to load CPWI data recorded on the Alpinion
%   Research scanner into a us_dataset class and beamform it with the USTB
%   routines. The data should be in the apropriate folder. The data can be
%   downloaded from: 
%   
%   https://www.dropbox.com/s/8sctpiyl6kir1w5/AlpinionCPWData.zip?dl=0
%
%   date : 17.02.2017
%   author: Ole Marius Hoel Rindal (olemarius@olemarius.net)

clear all
close all

% Set up filepath to files with parameters and data from the Alpinion scanner
%tag                 = '20161205-163945Angles29HyperCyst';
tag                 = '20161205-161341Angles29';

saveToHuff = 1;        %Do you want to save this dataset to Huff file?

parameter_filename  = ['data/Alpinion/CPW/',tag,'/Sequence.mat'];
data_filename       = ['data/Alpinion/CPW/',tag,'/40_layer0_idx40_BDATA_RF.mat']

%% initiate Alpinion object pointing to files with data
alp = alpinion();
alp.parameter_filename = parameter_filename;
alp.data_filename = data_filename

% initiate cpw_dataset object
cpw_dataset = cpw();
% read data from Alpinion 
alp.read(cpw_dataset);

%% convert to IQ data
cpw_dataset.demodulate(true,[],[],[],E.demodulation_algorithm.fieldsim);

% define a scan
scan=linear_scan();
scan.x_axis=linspace(cpw_dataset.geom(1,1),cpw_dataset.geom(end,1),256).';        % x vector [m]
scan.z_axis=linspace(0e-3,50e-3,256).';          % z vector [m]

F_number=1.75;
rx_beam_angle= 0;

orien=orientation();
orien.transmit_beam=beam(F_number,E.apodization_type.none);  
orien.receive_beam=beam(F_number,E.apodization_type.boxcar,rx_beam_angle);

%%
% define a reconstructions 
cpw_image=reconstruction();
cpw_image.scan=scan;
cpw_image.orientation=orien;

%% do beamforming
cpw_dataset.image_reconstruction(cpw_image,E.implementation.simple_matlab);
cpw_image.show();

if saveToHuff
    %% create a huff object
    recording=huff('alpinion_pw.h5','w');

    % adding datasets to staging area
    recording.stage(cpw_image);
    recording.stage(cpw_dataset); 

    
    recording.write();
end
