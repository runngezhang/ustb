%% Writting data to a UFF file
%
% In this example we show how to write channel and beamformed data into a
% UFF (Ultrasound File Format) file. The handling couldn't be simpler so
% this is going to be brief.
%
% _by Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no> 15.05.2017_


%% Getting channel data
%
% The first thing we need to save data into a UFF file is, you guessed it,
% data. Let us generate some channel data using the *fresnel*
% simulator included in the USTB. We won't get into details here. If you 
% want to know more about *fresnel* you can find some examples under the
% _fresnel_ folder.
%
% So here we define a 15 angles plane-wave sequence using a 128 elements
% linear array and a 5.2 MHz pulse. The phantom is a cross of point
% scatterers.

% phantom
x_sca=[zeros(1,7) -15e-3:5e-3:15e-3];
z_sca=[5e-3:5e-3:35e-3 20e-3*ones(1,7)];
N_sca=length(x_sca);
pha=uff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[x_sca.', zeros(N_sca,1), z_sca.', ones(N_sca,1)];    % point scatterer position [m]
             
% probe
prb=uff.linear_array();
prb.N=128;                  % number of elements 
prb.pitch=300e-6;           % probe pitch in azimuth [m]
prb.element_width=270e-6;   % element width [m]
prb.element_height=5000e-6; % element height [m]

% pulse
pul=uff.pulse();
pul.center_frequency=5.2e6;       % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]

% sequence
N_plane_waves=15; 
angles=linspace(-0.3,0.3,N_plane_waves);  % angle vector [rad]
seq=uff.wave();
for n=1:N_plane_waves 
    seq(n)=uff.wave();
    seq(n).probe=prb;
    seq(n).source.azimuth=angles(n);
    seq(n).source.distance=Inf;
    seq(n).sound_speed=pha.sound_speed;
end

% simulator
sim=fresnel();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% launch the simulation
channel_data=sim.go();

% setting dataset name & author information
channel_data.name = 'Test for UFF example';
channel_data.author = {'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Arun Nair <anair8@jhu.edu>'}; 
channel_data.reference = {'www.ustb.no'};

%% Getting beamformed data
%
% We will also generate some beamformed data to save into the same UFF
% file. To do that we define a scanning grid, a beamformer, and we set it
% to run.

% scan
scan=uff.linear_scan(linspace(-20e-3,20e-3,256).', linspace(0e-3,40e-3,256).');
 
% beamformer
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=scan;
bmf.receive_apodization.window=uff.window.tukey50;
bmf.receive_apodization.f_number=1.0;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uff.window.tukey50;
bmf.transmit_apodization.f_number=1.0;
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go({process.das_mex() process.coherent_compounding()});
b_data.plot();

%% Defining a UFF object
%
% Now that we have data to save we define a *uff* object to handle the UFF
% file. We do so by providing the full path (path + filename + extension)
% to the constructor of the *uff* class, for instance:

uff_file=uff('test02.uff');

%% 
%
% This will open 'test02.uff' file in the current folder. The constructor
% can take an additional parameter: _mode_ a string that specifies whether
% the file is meant to be for read-only ("read"), to read and write ("append"),
% or if we want to overwrite the file ("write"). By default the constructor
% opens the file in "append" mode.

%% Saving beamformed data
%
% It's about time we start saving some data. To do so we use the method
% *write* of the *uff* class. We must specify the object we want to save
% and the name it will have in the uff_data

uff_file.write(b_data,'b_data');

%%
%
% Now the beamformed data has been saves into the file. If you want to check
% the contents of the file with a HDF5 viewer such as
%
% <https://support.hdfgroup.org/products/java/release/download.html HDFView>
% 
% But we can check the contents of the file with the *index* method of
% *uff* with

display=true;
index=uff_file.index('/',display);

%% 
% 
% *uff/index* returns a cell with information on the datasets and
% groups in the specified location, see: 

index{:}

%%
% If the flag *display* is set then the
% function displays that information on screen. *uff/index* is not
% recursive: it only shows the contents of the specified location. 

%%
%
% If we try saving the data again with the same name...

try
    uff_file.write(b_data,'b_data');
catch me
    fprintf(2,[me.identifier ': ' me.message '\n']);
end

%%
%
% ...we get an error. Different datasets must have different names or be
% placed in different locations. For instance by:

uff_file.write(b_data,'b_data_copy');
uff_file.index('/',display);

%%
%
% It is also possible to save arrays of UFF structures. We can for instance
% define an array of beamformed data as

b_data_array=uff.beamformed_data();
for n=1:10
    b_data_array(n)=uff.beamformed_data();
    b_data_array(n).copy(b_data);
end

%%
%
% and store the whole array into the UFF file

uff_file.write(b_data_array,'b_data_array');
uff_file.index('/',display);

%% Saving channel data
% 
% Saving channel data (or any other *uff* structure) is exactly as with
% beamformed data. It might just take a bit more due to the larger amount
% of data. Here we save *scan* and *channel_data*

uff_file.write(scan,'scan');
uff_file.write(channel_data,'channel_data');
uff_file.index('/',display);
 
