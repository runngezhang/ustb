%% Computation of a CPWI dataset with Field II and beamforming with USTB
%
% This example shows how to load the data from a Field II simulation into a
% CPWI class, demodulate and beamform it with the USTB routines. The image
% is selected to have 3 different synthetic orientations, that are shown
% sequentially. The Field II simulation program (field-ii.dk) should be in 
% MATLAB's path.
%
% date:     11.03.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

%% basic constants
c0=1540;     % Speed of sound [m/s]
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

%% pulse definition
t0=(-1.0/bw/f0): dt : (1.0/bw/f0);
excitation=1;
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
noSubAz=round(width/(lambda/8));              % number of subelements in the azimuth direction
noSubEl=round(height/(lambda/8));             % number of subelements in the elevation direction
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
Na=5;                                       % number of plane waves 
alpha=linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]
F=10;                                       % number of frames

%% phantom
PRF=1./(2*40e-3/c0);             % pulse repetition frequency [Hz]
vz = c0*PRF/(8*f0*Na); vx=vz;   % velocity of the scatterer [m/s]
point=[0 0 20e-3]-[vx 0 -vz].*F/2.*Na/PRF;  % point initial position [m] 
cropat=round(2*2*40e-3/c0/dt);  % maximum time sample, samples after this will be dumped

%% output data
t_out=0:dt:((cropat-1)*dt);         % output time vector
CPW=zeros(cropat,N_elements,Na,F);  % impulse response channel data

%% Compute CPW signals
time_index=0;
disp('Field II: Computing CPW dataset');
wb = waitbar(0, 'Field II: Computing CPW dataset');
for f=1:F
    waitbar(0, wb, sprintf('Field II: Computing CPW dataset, frame %d',f));
    for n=1:Na
        waitbar(n/Na, wb);
        
        % position of the scatterer
        sca=point+[vx 0 -vz].*time_index/PRF;
        
        % transmit aperture
        xdc_apodization(Th,0,ones(1,N_elements));
        xdc_times_focus(Th,0,x0(:).'*sin(alpha(n))/c0);

        % receive aperture
        xdc_apodization(Rh, 0, ones(1,N_elements));
        xdc_focus_times(Rh, 0, zeros(1,N_elements));
        
        % do calculation
        [v,t]=calc_scat_multi(Th, Rh, sca, 1);

        % lag compensation
        t_in=(0:dt:((size(v,1)-1)*dt))+t-lag*dt;
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

%% Define a cpw dataset object
cpw_dataset=cpw('Field II, CPW, RF format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c0,...                            % reference speed of sound (m/s)
      alpha.',...                       % angle vector [rad]
      t_out.',...                       % time vector (s)
      PW,...                            % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,2)]);      % probe geometry [x, y, z] (m)

% demodulate signal
cpw_dataset.demodulate(true);

%% Define a scan object
scan=linear_scan();
scan.x_axis=linspace(-5.1e-3,5.1e-3,256).';              % x vector [m]
scan.z_axis=linspace(15e-3,25e-3,256).';                 % z vector [m]

%% Define a vector of orientations
F_number=1.75;
orientations=orientation();
orientations.transmit_beam=beam(1,E.apodization_type.none); % None apodization on transmit should be applied to produce low res images
orientations.receive_beam=beam(F_number,E.apodization_type.hamming,0);

%% Define a reconstruction object
recons=reconstruction();
recons.scan=scan;
recons.orientation=orientations;

%% Low res
raw_data=cpw_dataset.low_res(recons);
reshaped_data=reshape(raw_data,[size(recons.scan.x_matrix,1) size(recons.scan.x_matrix,2) size(cpw_dataset.data,3) size(cpw_dataset.data,4)]);

% show data
for f=1:size(reshaped_data,4)
    for n=1:size(reshaped_data,3)
        figure(101);
        imagesc(recons.scan.x_axis*1e3,recons.scan.z_axis*1e3,20*log10(abs(reshaped_data(:,:,n,f)))); colormap gray; axis equal; axis tight;
        xlabel('x[mm]'); ylabel('z[mm]');
        title(sprintf('angle=%d frame=%d',n,f));
        caxis([min(20*log10(abs(reshaped_data(:)))) max(20*log10(abs(reshaped_data(:))))])
        pause(0.25);
    end
end



