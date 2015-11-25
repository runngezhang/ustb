%% Comparing the different CPWI beamforming codes along time
% The file:
%
%       ps_cpw_iq.mat
%
% which must be located within MATLAB's path.
%
% date:     20.03.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

%% USTB 1.9

% Create a reconstruction object
recons=reconstruction();

% define the scan 
recons.scan=linear_scan();
recons.scan.x_axis=linspace(-20e-3,20e-3,256).';               % x vector [m]
recons.scan.z_axis=linspace(0e-3,40e-3,325).';                 % z vector [m]

% define the transmit & receive beams
F_number=1.1;
recons.orientation=orientation();
recons.orientation.transmit_beam=beam(F_number,E.apodization_type.none);
recons.orientation.receive_beam=beam(F_number,E.apodization_type.boxcar);

% format IQ
load('old_test_case.mat');                          

recons.name='CPW, IQ, Matlab';                                  % reconstruction name (optional)
tic
s.image_reconstruction(recons);  % request reconstruction
ustb1d9_time=toc;
recons.show();                                                % show 
title(sprintf('USTB 1.9 code: %0.4f ms',ustb1d9_time*1e3));

%% Thor Andreas
% p structure
p.c=s.c0;                                           % c         - Speed f sound [m/s]
p.dx=s.geom(2,2)-s.geom(2,1);                       % dx        - Array pitch
p.t0=s.time(1);                                     % t0        - Time of first sample
p.tx_angle=s.angle;                                 % tx_angle  - Transmit angle in radians
p.rx_angle=0.*s.angle;                              % rx_angle  - Receive angle in radians
p.FN=recons.orientation.receive_beam.f_number;      % FN        - F number for receive aperture.
p.fs_in=1./(s.time(2)-s.time(1));                    % fs_in     - Input Sampling frequency
p.fs_out=p.fs_in;                                   % fs_out    - Output Sampling frequency
p.f_demod=s.modulation_frequency;                   % f_demod   - IQ demodulation frequency. Needed if data is complex

% data
iq_ch_tx=s.data;
% scan data
x=recons.scan.x_axis;
r=recons.scan.z_axis.';
% apodization
apod2=single(reshape(s.receive_apodization,size(recons.scan.x_matrix,1),size(recons.scan.x_matrix,2),size(s.receive_apodization,2)));

% calling mex function
tic;
iq_bf0 = mex.planewave_beamforming2(single(iq_ch_tx),p,single(x),single(r),single(apod2));  % beamforming procedure
iq_bf0(isnan(iq_bf0)) = 0;
thor_time=toc;

recons.data=sum(iq_bf0,3);
recons.show(); 
title(sprintf('Thor Andreas code: %0.4f ms',thor_time*1e3));




