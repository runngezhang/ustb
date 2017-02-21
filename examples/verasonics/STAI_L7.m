%% Adquire and record a STA dataset with L7-4v probe

% date:     21.02.2017
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%           Ole Marius HOel Rindal <olemarius@olemarius.net>
%
% History:  Slight modification from the original STAI_L11.m script to fit
% the L7-4 probe we have at UiO.

% To run you should be in the Verasonics folder and activate it. For
% instance by:
%
% >> cd C:\Users\verasonics\Documents\Vantage-3.0.7
% >> activate
%
% The run and chose "Add to Path"
%
% To save the data:
%  
%  1.- Freeze
%  2.- Close the VSX window

clear all;
close all;

filename='a.mat';
folderdata=['data\' datestr(now,'yyyymmdd')];
mkdir(folderdata);            
filedata=['L7_STA_' datestr(now,'HHMMSS') '.h5'];
hufffile=[folderdata '\' filedata];
scan_area=[-19e-3 0e-3 19e-3 40e-3];
pixels=[255 255];

%% SI units
c0=1540;                % reference speed of sound [m/s]
f0=5.2e6;               % central frequency [Hz]
ex_cycles = 2.5;        % number of cycles of the excitation signal
ex_power = 0.67;        % signal duty cycle [0, 1] that relates to the amount of power delivered to the element  
ex_polarity = 1;        % easy way of changing the polarity
no_Frame = 1;           % number of frames to be acquired
no_Apert = 128;         % number of apertures per plane wave
PRF=9000;               % Pulse repetition frequency [pulses/s]
FN=1.2;                 % F-number
Fs=4*f0;                % sampling frequency

%% dependent values
lambda=c0/f0;                               % wavelength [m]

% Specify system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = c0;      % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.simulateMode = 1;       % 1 forces simulate mode, even if hardware is present.

%% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans.frequency = f0/1e6;   % The center frequency for the A/D 4xFc sampling.
% note nominal center frequency in computeTrans is 7.813 MHz
Trans = computeTrans(Trans);  % L12-3v transducer is 'known' transducer so we can use computeTrans.

%% Specify SFormat structure array.
SFormat.transducer = 'L7-4';   % 128 element linear array
SFormat.scanFormat = 'RLIN';     % rectangular linear array scan
SFormat.radius = 0;              % ROC for curved lin. or dist. to virt. apex
SFormat.theta = 0;
SFormat.numRays = 1;      % no. of Rays (1 for Flat Focus)
SFormat.FirstRayLoc = [0,0,0];   % x,y,z
SFormat.rayDelta = 128*Trans.spacing;  % spacing in radians(sector) or dist. between rays (wvlnghts)
SFormat.startDepth = 2;   % Acquisition depth in wavelengths
SFormat.endDepth = 3*128;%192;   % This should preferrably be a multiple of 128 samples.

%% Specify PData structure array.
PData.sFormat = 1;      % use first SFormat structure.
PData.Size(1) = pixels(2);  % Z
PData.Size(2) = pixels(1);  % X
PData.Size(3) = 1;      % single image page ???
PData.pdeltaX = ((scan_area(3)-scan_area(1))/PData.Size(2))/lambda;
PData.pdeltaZ = ((scan_area(4)-scan_area(2))/PData.Size(1))/lambda;
PData.Origin = [scan_area(1)/lambda, 0, scan_area(2)/lambda]; % x,y,z of upper lft crnr.

%% Specify Media object. 'pt1.m' script defines array of point targets.
Media.MP=[0, 0, 20e-3/lambda, 1; Trans.ElementPos];

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = no_Apert*4096; % this size allows for 3 acqs ???, maximum range
Resource.RcvBuffer(1).colsPerFrame = 256;
Resource.RcvBuffer(1).numFrames = no_Frame;     % 40 frames used for RF cineloop.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).rowsPerFrame = 2*1024;    % this is for greatest depth
Resource.InterBuffer(1).colsPerFrame = PData.Size(2);
Resource.InterBuffer(1).numFrames = 1;          % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).rowsPerFrame = 1024;
Resource.ImageBuffer(1).colsPerFrame = PData.Size(2);
Resource.ImageBuffer(1).numFrames = 10;         % 10 image buffer only
Resource.DisplayWindow(1).Title = 'L11-4v STAI';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).pdelta = 0.45;
Resource.DisplayWindow(1).Position = [250,250, ...    % upper left corner position
    ceil(PData.Size(2)*PData.pdeltaX/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData.Size(1)*PData.pdeltaZ/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),PData.Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Colormap = gray(256);

%% Specify Transmit waveform structure.  
TW.type = 'parametric';
% pulse expecification in [MHz, duty-cycle, number-of-half-cycles, boolean]
TW.Parameters = [f0/1e6, ex_power, ex_cycles*2, ex_polarity];

%% Specify TX structure array.  
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, no_Apert); % number of apertures + 2 extra dummy transmits
% writing sequence angles
for n = 1:no_Apert
    TX(n).Apod(n) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n)); 
end
              
%% Specify TGC Waveform structure.
TGC.CntrlPts = [139,535,650,710,770,932,992,1012];
TGC.rangeMax = SFormat.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
maxAcqLength = sqrt(SFormat.endDepth^2 + (Trans.numelements*Trans.spacing)^2) - SFormat.startDepth;
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,128), ...
                        'startDepth', SFormat.startDepth, ...
                        'endDepth', SFormat.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'samplesPerWave', 4, ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,no_Apert*no_Frame);
%% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames  % no_Apert*acquisitions per frame
    k = no_Apert*(i-1);
    for n = 1:no_Apert; % for each aperture acquire all angles    
        Receive(k+n).framenum = i;
        Receive(k+n).acqNum = n;
    end
end

%% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', (1:no_Apert)');
% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 4, ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, no_Apert);
% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 3;
for n=1:no_Apert
    ReconInfo(n).rcvnum = n;
    ReconInfo(n).txnum = n;  
end

ReconInfo(no_Apert).mode = 5;

%% Specify Process structure array.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',2.0,...            % pgain is image processing gain
                         'persistMethod','simple',...
                         'persistLevel',30,...
                         'interp',1,...      % method of interpolation (1=4pt interp)
                         'compression',0.5,...      % X^0.5 normalized to output word size
                         'reject',2,...
                         'mappingMode','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

%% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 1/PRF*1e6;       % PRF
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 1/PRF*1e6; % framerate
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'timeToNextAcq';  % time between ???
SeqControl(5).argument = 10;              % ???
nsc = 6; % nsc is count of SeqControl objects

%% Specify Event structure arrays.
n = 1; % n is count of Events
for i = 1:Resource.RcvBuffer(1).numFrames
    k = no_Apert*(i-1);
    for j = 1:no_Apert
        % first transmit - first receive
        Event(n).tx = j;         % use 1st TX structure.
        Event(n).rcv = k+j;      % use 1st Rcv structure.
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2; % time between syn. aper. acqs.
        n = n+1;
    end
    % Replace last Event's seqControl value.
    Event(n-1).seqControl = [3,nsc]; % time between frames, SeqControl struct defined below.
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    
    Event(n).info = 'Reconstruct & process'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % processing
    Event(n).seqControl = 0;
    if floor(i/2) == i/2     % Exit to Matlab every 4th frame 
        Event(n).seqControl = 4;
    end
    n = n+1;
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 1; % jump command


%% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%-UI#1Callback');

% - Range Change
UI(2).Control = {'UserA1','Style','VsSlider','Label','Range',...
                 'SliderMinMaxVal',[64,320,SFormat.endDepth],'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%-UI#2Callback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = no_Apert;

% Save all the structures to a .mat file.
save('a');

%% call VSX
VSX;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% converting the format to USTB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Converting data format to USTB');
no_samples=Receive(1).endSample;
data=zeros(no_samples, 128, 128, Resource.RcvBuffer(1).numFrames);
geom=Trans.ElementPos(1:128,1:3)*1e-3;     
delay_x0=sqrt(sum((geom-ones(128,1)*[0 0 20e-3]).^2,2))/c0;

% offset calculation
offset_distance=(TW.peak)*lambda;   % in [m]
if strcmp(Trans.units,'mm')
    offset_distance=offset_distance+2*Trans.lensCorrection*1e-3;
elseif strcmp(Trans.units,'wavelengths')
    offset_distance=offset_distance+2*Trans.lensCorrection*lambda;
end
offset_time=offset_distance/Resource.Parameters.speedOfSound;   % in [s]

% convert data
n=1;
t_out=0:(1/Fs):((no_samples-1)/Fs);
plot_delayed_signal=0;
for n_frame = 1:Resource.RcvBuffer(1).numFrames
    for n_tx = 1:no_Apert
        % compute time vector for this line
        t_ini=2*Receive(n).startDepth*lambda/c0;
        t_end=2*Receive(n).endDepth*lambda/c0;
        no_t=(Receive(n).endSample-Receive(n).startSample+1);
        t_in=linspace(t_ini,t_end,no_t)-offset_time;%scanOffsetTime;
        
        % read data
        data(:,:,n_tx,n_frame)=interp1(t_in,double(RcvData{1}(Receive(n).startSample:Receive(n).endSample,:,n_frame)),t_out,'linear',0);
        n=n+1;

        % to check delay calculation
        if plot_delayed_signal
            delay= delay_x0+delay_x0(n_tx);

            figure(101); hold off;
            pcolor(1:no_Apert,t_out,real(data(:,:,n_tx,n_frame))); shading flat; colormap gray; colorbar; hold on;
            plot(1:no_Apert,delay,'r');
            title(n_tx);
            ylim([0.9*min(delay) 1.1*max(delay)]);
            pause();
        end
    end
end

%% create USTB data class structure
sta_dataset=sta('STA Verasonics',E.signal_format.RF,c0,t_out.',double(data),geom);

%% save RF as huff
disp('Writting HUFF file');
huff_file=huff(hufffile,'w');
huff_file.stage(sta_dataset); 
huff_file.write();

%% demodulate data
sta_dataset.demodulate(true,[],[0 1e6 1.9*f0 2*f0],[],E.demodulation_algorithm.fastfon);

%% image reconstruction with USTB 
sta_image=reconstruction();

% define the scan 
sta_image.scan=linear_scan();
sta_image.scan.x_axis=linspace(-5e-3,5e-3,256).';                 % x vector [m]
sta_image.scan.z_axis=linspace(15e-3,25e-3,256).';                  % z vector [m]

% define the transmit & receive beams
sta_image.orientation=orientation();
sta_image.orientation.transmit_beam=beam(1.7,E.apodization_type.boxcar);
sta_image.orientation.receive_beam=beam(1.2,E.apodization_type.boxcar);

sta_image.name='STA L11 example';                                   % reconstruction name (optional)
sta_dataset.image_reconstruction(sta_image);                        % request reconstruction
im_1=sta_image.show('log',60);                                      % show 
    
return


% **** Callback routines to be converted by text2cell function. ****
%-UI#1Callback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%-UI#1Callback

%-UI#2Callback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','SFormat.endDepth'));
    return
end
range = UIValue;
assignin('base','range',range);
SFormat = evalin('base','SFormat');
SFormat.endDepth = range;
assignin('base','SFormat',SFormat);
evalin('base','PData.Size(1) = ceil((SFormat.endDepth-SFormat.startDepth)/PData.pdeltaZ);');
evalin('base','[PData.Region,PData.numRegions] = createRegions(PData);');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.pdeltaZ/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
Trans = evalin('base', 'Trans');
maxAcqLength = sqrt(range^2 + (Trans.numelements*Trans.spacing)^2)-SFormat.startDepth;
wlsPer128 = 128/(4*2);
for i = 1:size(Receive,2)
    Receive(i).endDepth = SFormat.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128);
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = SFormat.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'SFormat','PData','Receive','Recon','DisplayWindow','ImageBuffer'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%-UI#2Callback
