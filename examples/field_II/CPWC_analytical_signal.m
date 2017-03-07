%% Demonstration of "the analytical signal" on a CPWI dataset from a Field II simulation
%
% This example shows how "the analytical signal (AS)" can be used. The
% analytical signal is a complex version of the RF signal with the original
% RF signal as the real part, and the Hilbert transform as the imaginary
% part. 
%
% NB! The Matlab function Hilbert() returns the analytical signal
%
% See appendix B in http://www.olemarius.net/Thesis/thesis.pdf for more
% info.
%
% The data is created with a Field II simulation saved into CPWI class, 
% one using RF one using AS and one using IQ. We show that by using AS we
% avoid the problem of having a oversampled image in the z-direction.
%
% The Field II simulation program (field-ii.dk) should be in MATLAB's path.
%
% date:     06.03.2017
% authors:  Ole Marius Hoel Rindal <olemarius@olemarius.net>
%           Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>      

close all;
clear all;

%% basic constants
c0=1540;    % Speed of sound [m/s]
fs=100e6;   % Sampling frequency [Hz]
dt=1/fs;    % Sampling step [s] 

%% field II initialisation
field_init;
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% transducer definition L9-4/38 Ultrasonix
f0=5e6;             % Transducer center frequency [Hz]
bw=0.65;            % probe bandwidth [1]
lambda=c0/f0;       % Wavelength [m]
height=6e-3;        % Height of element [m]
pitch=0.3048e-3;    % pitch [m]
kerf=0.035e-3;      % gap between elements [m]
width=pitch-kerf;   % Width of element [m]
lens_el=19e-3;      % position of the elevation focus
N_elements=128;     % Number of elements
pulse_duration=2.5; % pulse duration

%% pulse definition
t0=(-1.0/bw/f0): dt : (1.0/bw/f0);
te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation=square(2*pi*f0*te+pi/2);
impulse_response=gauspuls(t0, f0, bw);
two_ways_ir= conv(conv(impulse_response,impulse_response),excitation);
if mod(length(impulse_response),2)
    lag=(length(two_ways_ir)-1)/2;          
else
    lag=(length(two_ways_ir))/2;
end

% show the pulse to check that the lag estimation is on place (and that the pulse is symmetric)
figure;
plot((0:(length(two_ways_ir)-1))*dt -lag*dt,two_ways_ir); hold on; grid on; axis tight
plot((0:(length(two_ways_ir)-1))*dt -lag*dt,abs(hilbert(two_ways_ir)),'r')
plot([0 0],[min(two_ways_ir) max(two_ways_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');

%% aperture objects
% definition of the mesh geometry
noSubAz=2;              % number of subelements in the azimuth direction
noSubEl=8;              % number of subelements in the elevation direction
Th = xdc_linear_array (N_elements, width, height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
Rh = xdc_linear_array (N_elements, width, height, kerf, noSubAz, noSubEl, [0 0 Inf]); 

% setting excitation, impulse response and baffle
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);

%% get geometrical center of elements
data = xdc_get(Th,'rect');
geo=data(24:26,:);
x0=zeros(1,N_elements);
z0=zeros(1,N_elements);
for n=1:N_elements
    n_ini=(n-1)*noSubAz*noSubEl+1;
    n_fin=n_ini+noSubAz*noSubEl-1;
    x0(n)=mean(geo(1,n_ini:n_fin));
    z0(n)=mean(geo(3,n_ini:n_fin));
end

%% plane wave sequence
alpha_max=15*pi/180;                        % maximum angle [rad]
Na=1;                                      % number of plane waves 
if (Na>1)
    alpha=linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]
else
    alpha=0;
end
F=1;                                        % number of frames

%% phantom
PRF=1./(2*40e-3/c0);                                % pulse repetition frequency [Hz]
x_point=[zeros(1,8) linspace(-20e-3,20e-3,9)];                  
z_point=[linspace(5e-3,40e-3,8) 20e-3*ones(1,9)];

point=[x_point.' zeros(length(x_point),1) z_point.'];            % point initial position [m] 
cropat=round(8*max(sqrt(sum(point.^2,2)))/c0/dt);   % maximum time sample, samples after this will be dumped

%% output data
t_out=0:dt:((cropat-1)*dt);         % output time vector
CPW=zeros(cropat,N_elements,Na,F);  % impulse response channel data

%% Compute CPW signals
sca=point;
time_index=0;
disp('Field II: Computing CPW dataset');
wb = waitbar(0, 'Field II: Computing CPW dataset');
for f=1:F
    waitbar(0, wb, sprintf('Field II: Computing CPW dataset, frame %d',f));
    for n=1:Na
        waitbar(n/Na, wb);
        
        % transmit aperture
        xdc_apodization(Th,0,ones(1,N_elements));
        xdc_times_focus(Th,0,x0(:).'*sin(alpha(n))/c0);

        % receive aperture
        xdc_apodization(Rh, 0, ones(1,N_elements));
        xdc_focus_times(Rh, 0, zeros(1,N_elements));
        
        % do calculation
        [v,t]=calc_scat_multi(Th, Rh, sca, ones(size(point,1),1));

        % lag compensation
        t_in=(0:dt:((size(v,1)-1)*dt))+t;%-lag*dt;
        v_aux=interp1(t_in,v,t_out,'linear',0);

        % build the dataset
        CPW(:,:,n,f)=v_aux;
        
        time_index=time_index+1;
    end
end
close(wb);

%% USTB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define a cpw dataset object for the IQ data
cpw_dataset_iq=cpw('Field II, CPW, IQ format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c0,...                            % reference speed of sound (m/s)
      alpha.',...                       % angle vector [rad]
      t_out.',...                       % time vector (s)
      CPW*10^24,...                           % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,1) z0.']);      % probe geometry [x, y, z] (m)

% convert to IQ data
cpw_dataset_iq.demodulate(true,[],[],[],E.demodulation_algorithm.fieldsim);  
  
% Define a cpw dataset object for the AS data
cpw_dataset_as=cpw('Field II, CPW, AS format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c0,...                            % reference speed of sound (m/s)
      alpha.',...                       % angle vector [rad]
      t_out.',...                       % time vector (s)
      CPW*10^24,...                           % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,1) z0.']);      % probe geometry [x, y, z] (m)

% convert to AS data
cpw_dataset_as.make_analytical_signal();

% Define a cpw dataset object for the RF data
cpw_dataset_rf=cpw('Field II, CPW, RF format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c0,...                            % reference speed of sound (m/s)
      alpha.',...                       % angle vector [rad]
      t_out.',...                       % time vector (s)
      CPW*10^24,...                           % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,1) z0.']);      % probe geometry [x, y, z] (m)



%% define a scan
scan=linear_scan();
scan.x_axis=linspace(-20e-3,20e-3,256).';        % x vector [m]
scan.z_axis=linspace(0e-3,40e-3,256).';          % z vector [m]

% define a synthetic orientation
F_number=1.75;
rx_beam_angle=[0]*pi/180;

orien(1)=orientation();
orien(1).transmit_beam=beam(F_number,E.apodization_type.none);  
orien(1).receive_beam=beam(F_number,E.apodization_type.boxcar,rx_beam_angle(n));


% define a reconstruction for each datatype
cpw_image_iq=reconstruction();
cpw_image_iq.scan=scan;
cpw_image_iq.orientation=orien;

cpw_image_as=reconstruction();
cpw_image_as.scan=scan;
cpw_image_as.orientation=orien;

cpw_image_rf=reconstruction();
cpw_image_rf.scan=scan;
cpw_image_rf.orientation=orien;

%% beamform and show images
cpw_dataset_iq.image_reconstruction(cpw_image_iq,E.implementation.simple_matlab);
cpw_image_iq.show('log',60,'Beamformed from RF data');

cpw_dataset_as.image_reconstruction(cpw_image_as,E.implementation.simple_matlab);
cpw_image_as.show('log',60,'Beamformed from AS data');

cpw_dataset_rf.image_reconstruction(cpw_image_rf,E.implementation.simple_matlab);
cpw_image_rf.show('log',60,'Beamformed from RF data');    

% We demonstrate that both AS and IQ have no aliasing in the
% envelope detection and thus no aliasing in the image. While the image
% created from RF gets aliasing because of an inaccurate Hilbert transform
% when detecting the envelope. Therefore it might be wise to use AS instead
% of RF data.

%% Plot the Fourier Transform of the different data to show the difference

figure
subplot(311)
rf_x_axis = linspace(-cpw_dataset_rf.sampling_frequency/2,...
             cpw_dataset_rf.sampling_frequency/2,...
             length(cpw_dataset_rf.data(:,end/2)));
plot(rf_x_axis*10^-6,db(abs(fftshift(fft(cpw_dataset_rf.data(:,end/2))))))
ax(1) = gca; xlabel('Frequency [MHz]'); ylabel('Amplitude [dB]');
title('FFT of RF signal');

subplot(312)
as_x_axis = linspace(-cpw_dataset_as.sampling_frequency/2,...
             cpw_dataset_as.sampling_frequency/2,...
             length(cpw_dataset_as.data(:,end/2)));
plot(as_x_axis*10^-6,db(abs(fftshift(fft(cpw_dataset_as.data(:,end/2))))))
ax(2) = gca; xlabel('Frequency [MHz]'); ylabel('Amplitude [dB]');
title('FFT of AS signal');

subplot(313)
iq_x_axis = linspace(-cpw_dataset_iq.sampling_frequency/2,...
             cpw_dataset_iq.sampling_frequency/2,...
             length(cpw_dataset_iq.data(:,end/2)));
plot(iq_x_axis*10^-6,db(abs(fftshift(fft(cpw_dataset_iq.data(:,end/2))))))
ax(3) = gca; xlabel('Frequency [MHz]'); ylabel('Amplitude [dB]');
title('FFT of IQ signal');


linkaxes(ax);
ylim([-200 100]);
xlim([-cpw_dataset_rf.sampling_frequency*10^-6/4 cpw_dataset_rf.sampling_frequency*10^-6/4]);
  
%% Plot of data from one channel for each type of data
figure
subplot(311)
x_axis = 1:length(cpw_dataset_rf.data(:,end/2,:,:,1));
plot(x_axis,cpw_dataset_rf.data(:,end/2,:,:,1));
ax(1) = gca;
xlabel('Samples');
ylabel('Amplitude');
title('One channel of RF data');

subplot(312)
plot(x_axis,real(cpw_dataset_as.data(:,end/2,:,:,1)),'DisplayName','real');hold on
plot(x_axis,imag(cpw_dataset_as.data(:,end/2,:,:,1)),'DisplayName','imag','Color','r');
ax(2) = gca;
legend show
xlabel('Samples');
ylabel('Amplitude');
title('One channel of AS data');

subplot(313)
x_axis_iq = linspace(x_axis(1),x_axis(end),length(cpw_dataset_iq.data(:,end/2,:,:,1)));
plot(x_axis_iq,real(cpw_dataset_iq.data(:,end/2,:,:,1)),'DisplayName','real');hold on
plot(x_axis_iq,imag(cpw_dataset_iq.data(:,end/2,:,:,1)),'DisplayName','imag','Color','r');
ax(3) = gca;
xlabel('Samples NB! Be aware that IQ have a factor four less samples');
ylabel('Amplitude');
title('One channel of IQ data');
legend show
linkaxes(ax);xlim([1800 3200]);



%% Plot of data from one depth line (z-direction) in the final image for each type of data
figure
subplot(311)
x_axis = 1:length(cpw_image_rf.data(:,end/2,:,:,1));
plot(x_axis,cpw_image_rf.data(:,end/2,:,:,1));
ax(1) = gca;
xlabel('Samples');
ylabel('Amplitude');
title('One axial line in the image of RF data');

subplot(312)
plot(x_axis,real(cpw_image_as.data(:,end/2,:,:,1)),'DisplayName','real');hold on
plot(x_axis,imag(cpw_image_as.data(:,end/2,:,:,1)),'DisplayName','imag','Color','r');
ax(2) = gca;
legend show
xlabel('Samples');
ylabel('Amplitude');
title('One axial line in the image of AS data');

subplot(313)
x_axis_iq = linspace(x_axis(1),x_axis(end),length(cpw_image_iq.data(:,end/2,:,:,1)));
plot(x_axis_iq,real(cpw_image_iq.data(:,end/2,:,:,1)),'DisplayName','real');hold on
plot(x_axis_iq,imag(cpw_image_iq.data(:,end/2,:,:,1)),'DisplayName','imag','Color','r');
ax(3) = gca;
xlabel('Samples NB! Be aware that IQ have a factor four less samples');
ylabel('Amplitude');
title('One axial line in the image of IQ data');
legend show
linkaxes(ax);%xlim([1800 3200]);

