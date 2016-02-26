%% Speedtest of CPWC mex implementations

% date:     11.03.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

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
for n=1:N_elements
    n_ini=(n-1)*noSubAz*noSubEl+1;
    n_fin=n_ini+noSubAz*noSubEl-1;
    x0(n)=mean(geo(1,n_ini:n_fin));
end

%% plane wave sequence
alpha_max=15*pi/180;                        % maximum angle [rad]
Na=3;                                      % number of plane waves 
if (Na>1)
    alpha=linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]
else
    alpha=0;
end
F=2;                                        % number of frames

%% phantom
PRF=1./(2*40e-3/c0);                                % pulse repetition frequency [Hz]
x_point=[zeros(1,8) linspace(-20e-3,20e-3,9)];                  
z_point=[linspace(5e-3,40e-3,8) 20e-3*ones(1,9)];
%[X Z]=meshgrid(x_point,z_point);
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

%CPW=repmat(CPW,[1 1 1 100]);

%% USTB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define a cpw dataset object
cpw_dataset=cpw('Field II, CPW, RF format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c0,...                            % reference speed of sound (m/s)
      alpha.',...                       % angle vector [rad]
      t_out.',...                       % time vector (s)
      CPW,...                           % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,2)]);      % probe geometry [x, y, z] (m)

% convert to IQ data
cpw_dataset.demodulate(true,[],[],[],E.demodulation_algorithm.fieldsim);

% define a scan
scan=linear_scan();
scan.x_axis=linspace(-20e-3,20e-3,256).';        % x vector [m]
scan.z_axis=linspace(0e-3,40e-3,256).';        % z vector [m]

% define a synthetic orientation
F_number=1.75;
orien(1)=orientation();
orien(1).transmit_beam=beam(F_number,E.apodization_type.none);  
orien(1).receive_beam=beam(F_number,E.apodization_type.boxcar,0*pi/180);
orien(2)=orientation();
orien(2).transmit_beam=beam(F_number,E.apodization_type.none);  
orien(2).receive_beam=beam(F_number,E.apodization_type.boxcar,-10*pi/180);
orien(3)=orientation();
orien(3).transmit_beam=beam(F_number,E.apodization_type.none);  
orien(3).receive_beam=beam(F_number,E.apodization_type.boxcar,10*pi/180);

% define a reconstructions 
cpw_image_ustb_mex=reconstruction();
cpw_image_ustb_mex.scan=scan;
cpw_image_ustb_mex.orientation=orien;

cpw_image_ta=reconstruction();
cpw_image_ta.scan=scan;
cpw_image_ta.orientation=orien;

cpw_image_lr=reconstruction();
cpw_image_lr.scan=scan;
cpw_image_lr.orientation=orien;

% convert to IQ data
cpw_dataset.data=repmat(cpw_dataset.data(:,:,:,1),[1 1 1 f]);
    
%% beamform ustb mex
mex_time(f)=cpw_dataset.image_reconstruction(cpw_image_ustb_mex,E.implementation.mex);
cpw_image_ustb_mex.show();

%% beamform thor-andreas code
ta_time(f)=cpw_dataset.image_reconstruction(cpw_image_ta,E.implementation.thor_andreas);
cpw_image_ta.show();

%% beamform low-resolution code
lr_time(f)=cpw_dataset.image_reconstruction(cpw_image_lr,E.implementation.low_resolution);
cpw_image_lr.show();
