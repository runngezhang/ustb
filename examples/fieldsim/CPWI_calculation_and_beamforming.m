%% Simulation with FieldSim v3 and beamforming with USTB
%
% Example on how to launch a Field II simulation with FieldSim and
% beamform the data with USTB. The FielsSim toolbox must be installed 
% and the Field II program must be included into MATLAB's path. To run this 
% example the following data file is needed
%
%       PhantomSimple.m
%
% which must be located in the corresponding FieldSim folder.
%
% date:     11.03.2015
% authors:  Ingvild Kinn Ekroll <ingvild.k.ekroll@ntnu.no>
%           Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

%% Start FieldSim
sim           = FieldSim.Simulation();
sim.name      = 'PSF-test-beamforming';
sim.simulator = 'Field2ChannelData';
doIQdemodulation = 0; % '1', '0': %If sim.simulator is Field2ChannelData (giving channel data), '1' keeps data in RF domain before beamforming --> gives beamformed RF data

%% Select the FieldSim basic configuration
sim.probe = '9L';
sim.selectMode('BmodeChannelData');
sim.scan = 'BmodePlaneWave';

%% Configure Simulation
% specification of the transmitted pulse 
sim.fs = 100e6;                                 % sampling frequency [Hz]
sim.scan.txPulser.f0 = 6e6;                     % pulse frequency [Hz]
sim.scan.txPulser.noPeriodsExcitation = 2.5;    % pulse duration [cycles]

% specification of the plane-wave sequence
sim.scan.maxAngle = 30*pi/180;                 % angle span [rad]  CAUTION: This is not the maximum angle but the angle span covered
sim.scan.noTxBeams = 5;                        % number of plane waves
sim.scan.txApod.dynamicAperture=0;
sim.scan.txApod.apodization={'none','none'};
sim.scan.rxApod.dynamicAperture=0;
sim.scan.rxApod.apodization={'none','none'};

% specification of the scanned region
sim.scan.scanShape.range_min = 0e-3;  % minimum depth [m]
sim.scan.scanShape.range_max = 40e-3; % maximum depth [m]

% plot scan configuration
figure(), clf, hold all, grid on
plot(sim.scan);
axis equal tight
elc = 1e3*sim.probe.getElementCenters();
elc = squeeze(elc(:,round(end/2),:));
plot(elc(:,1), elc(:,3), 'bx')
ax_x = get(gca, 'XLim'); ax_z = get(gca, 'ZLim');

%% Doppler specifics
sim.scan.noPackets = 1;                         % number of packects
sim.scan.noFrames = 1;                          % number of frames
sim.scan.PRF = 4000*sim.scan.noTxBeams(1);      % PRF [Hz]

%% Phantom definition
z0=20e-3;                          % depth of the point scatterer [m]
sim.phantom = 'PhantomSimple';     % selection of the FieldSim phantom
sim.phantom.noScatterers=1;        % number of scatterers in a line
sim.phantom.rStart= z0;            % start depth for the line of scatterers [m]
sim.phantom.rEnd= z0;              % end depth for the line of scatterers [m]

%% Lag calculation 
excitation=sim.scan.txPulser.generateExcitation(sim.fs);    % this is to get the excitation length
ir = sim.probe.impulseResponse.getResampled(sim.fs);        % this is to get the impulse response length
two_ways_ir= conv(conv(ir,ir),excitation);                  % two ways pulse
lag=length(two_ways_ir)/2;                                  % the lag equals the sample in the middle

% show the pulse, just to check that the estimation is on place (and that the pulse is symmetric)
figure;
plot((0:(length(two_ways_ir)-1))/sim.fs -lag/sim.fs,two_ways_ir); hold on; grid on; axis tight
plot((0:(length(two_ways_ir)-1))/sim.fs -lag/sim.fs,abs(hilbert(two_ways_ir)),'r')
plot([0 0],[min(two_ways_ir) max(two_ways_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response FieldSim');

sim.simulator.autoFilterLagCorrection=0;         % Activate/Deactivate the lag correction in FieldSim. 
sim.simulator.manualFilterLagOffset=-lag/sim.fs; % insert our lag correction

%% Do simulation
sim.simulator.doScan();                          % do the simulation

%% converting data to USTB format

% time vector
time=sim.scan.startTime+(0:(size(sim.data.channel_data,1)-1)).'/sim.fs;

% angle vector
angles=zeros(length(sim.scan.scanSeq),1);
for kk = 1:length(sim.scan.scanSeq)
    angles(kk)=sim.scan.scanSeq{kk}.txBeam.tilt(1,1);
end

% data
data=reshape(sim.data.channel_data,[size(sim.data.channel_data,1) size(sim.data.channel_data,2) size(sim.data.channel_data,4) size(sim.data.channel_data,7)]);

%% define the dataset and demodulate
cpw_dataset=cpw('Ingvild data, CPW',...                   % dataset name
                 E.signal_format.RF,...                   % signal format (RF/IQ)
                 sim.propagation.c,...                    % reference speed of sound (m/s)
                 angles,...                               % vector of plane wave angles (rad)
                 time,...                                 % time vector (s)
                 data,...                                 % data [samples,channels,firings,frames]
                 squeeze(sim.probe.getElementCenters())); % probe geometry [x y z] (m)

% % show RF data
% figure(200);
% for n=1:length(cpw_dataset.angle)
%     imagesc(1:size(cpw_dataset.data,2),cpw_dataset.time,abs(cpw_dataset.data(:,:,n))); colormap gray; hold on;
%     delay=(z0*cos(cpw_dataset.angle(n))+sqrt(cpw_dataset.geom(:,1).^2+z0^2))/cpw_dataset.c0;
%     plot(1:size(cpw_dataset.data,2),delay,'r--');
%     axis([1 size(cpw_dataset.data,2) min(delay) max(delay)]);
%     title(cpw_dataset.angle(n));
%     pause();
% end
             
% demodulate
cpw_dataset.demodulate(false,6e6);   

%% Define a reconstruction object
recons_fieldsim=reconstruction();

% define a linear scan 
recons_fieldsim.scan=linear_scan();
recons_fieldsim.scan.x_axis=linspace(-1.5e-2,1.5e-2,512).';               % x vector [m]
recons_fieldsim.scan.z_axis=linspace(15e-3,25e-3,512).';              % z vector [m]

% define the transmit & receive beams
%F-number, transmit apodization, steering angle [rad], length of the edge smoothing area [elements]
recons_fieldsim.orientation=orientation();
recons_fieldsim.orientation.transmit_beam=beam(2,E.apodization_type.boxcar,0,0); 
recons_fieldsim.orientation.receive_beam=beam(2,E.apodization_type.boxcar,0,0);

% request reconstruction 
cpw_dataset.image_reconstruction(recons_fieldsim);

% show it
recons_fieldsim.show();