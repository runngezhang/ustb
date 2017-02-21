%% Reading Verasonics data with HUFF 
%
% This example shows how you can read Verasonics data using HUFF 

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 12/02/2017$

%% data location
url='http://hirse.medisin.ntnu.no/ustb/data/verasonics/';   % if not found data will be downloaded from here
local_path='data/verasonics/';                              % location of example data in this computer                      

%% Create a reconstruction object
r=reconstruction();

% define the scan 
r.scan=linear_scan();
r.scan.x_axis=linspace(-5e-3,5e-3,256).';                  % x vector [m]
r.scan.z_axis=linspace(15e-3,25e-3,256).';                 % z vector [m]

% define the transmit & receive beams
r.orientation=orientation();
r.orientation.transmit_beam=beam(1.7,E.apodization_type.boxcar);
r.orientation.receive_beam=beam(1.2,E.apodization_type.boxcar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STA Verasonics simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if data is available & download
filename='L11_STA_simulation.h5';
a=dir([local_path filename]);
if not(numel(a))
    disp(['Downloading example data from ' url '. This may take a bit.']);
    mkdir(local_path);
    urlwrite([url filename],[local_path filename]);
end

% create a huff object
recording=huff([local_path filename],'r');

% open file and load dataset
recording.read();
s=recording.ultrasound_signal{1};

% demodulate
s.demodulate();

% reconstruct 
s.image_reconstruction(r);
r.show();
title('STA simulation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STA Verasonics experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if data is available & download
filename='L11_STA_phantom.h5';
a=dir([local_path filename]);
if not(numel(a))
    disp(['Downloading example data from ' url '. This may take a bit.']);
    mkdir(local_path);
    urlwrite([url filename],[local_path filename]);
end

% create a huff object
recording=huff([local_path filename],'r');

% open file and load dataset
recording.read();
s=recording.ultrasound_signal{1};

% demodulate
s.demodulate();

% reconstruct 
s.image_reconstruction(r);
r.show();
title('STA phantom');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CPW Verasonics simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if data is available & download
filename='L11_CPW_simulation.h5';
a=dir([local_path filename]);
if not(numel(a))
    disp(['Downloading example data from ' url '. This may take a bit.']);
    mkdir(local_path);
    urlwrite([url filename],[local_path filename]);
end

% create a huff object
recording=huff([local_path filename],'r');

% open file and load dataset
recording.read();
s=recording.ultrasound_signal{1};

% demodulate
s.demodulate();

% reconstruct 
s.image_reconstruction(r);
r.show();
title('CPWC simulation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CPW Verasonics experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if data is available & download
filename='L11_CPW_phantom.h5';
a=dir([local_path filename]);
if not(numel(a))
    disp(['Downloading example data from ' url '. This may take a bit.']);
    mkdir(local_path);
    urlwrite([url filename],[local_path filename]);
end

% create a huff object
recording=huff([local_path filename],'r');

% open file and load dataset
recording.read();
s=recording.ultrasound_signal{1};

% demodulate
s.demodulate();

% reconstruct 
s.image_reconstruction(r);
r.show();
title('CPWC phantom');