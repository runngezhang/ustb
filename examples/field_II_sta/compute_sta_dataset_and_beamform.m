clear all;
close all;

%% field II initialisation
field_init;
set_field('att', 0);
c=1540; set_field('c',c);                % Speed of sound [m/s]
fs=100e6; set_field('fs',fs); dt=1/fs;    % Sampling frequency [Hz]
set_field('use_rectangles',1);

% define transducer
f0=4.5e6; % Transducer center frequency [Hz]
lambda=c/f0; % Wavelength [m]
height=4e-3;  % Height of element [m]
pitch=0.33e-3; % pitch [m]
kerf=pitch./10; % Distance between transducer elements [m]
width=pitch-kerf; % Width of element [m]
N_elements=128; % Number of elements

% defining transmit and receive apertures
Ny=round(height/(lambda/8));
Nx=round(width/(lambda/8));
Th = xdc_linear_array (N_elements, width, height, kerf, Nx, Ny, [0 0 10000]);
Rh = xdc_linear_array (N_elements, width, height, kerf, Nx, Ny, [0 0 10000]);
data = xdc_get(Th,'rect');
geo=data(24:26,:);
x0=zeros(1,N_elements);
for n=1:N_elements
    n_ini=(n-1)*Nx*Ny+1;
    n_fin=n_ini+Nx*Ny-1;
    x0(n)=mean(geo(1,n_ini:n_fin));
end

% setting minimal delay 
min_delay=1e-7;
xdc_focus_times(Th, 0, min_delay*ones(1,N_elements));
xdc_focus_times(Rh, 0, min_delay*ones(1,N_elements));

% define excitation
cycles=5;
t0=0:dt:(cycles/f0); t0=t0-cycles/f0/2;
excitation=[1 zeros(1,length(t0))];
impulse_response=sin(2*pi*f0*t0).*exp(-0.5*t0.^2/(1/f0/2.25)^2); % gaussian pulse
xdc_impulse (Th, impulse_response);
xdc_excitation (Th, excitation);
xdc_impulse (Rh, impulse_response);

% define two-ways pulse
ir_2ways=conv(impulse_response,impulse_response); ir_2ways=ir_2ways./max(abs(ir_2ways));
t0_2ways=dt:dt:(dt*(length(ir_2ways))); t0_2ways=t0_2ways-2*cycles/f0/2-dt;

% Define domain
x_point=-15e-3:5e-3:15e-3;
z_point=5e-3:5e-3:40e-3;
[xxp zzp]=meshgrid(x_point,z_point);
N_sca=length(zzp(:));
point=[xxp(:) zeros(N_sca,1) zzp(:)];

cropat=round(2.5*2*max(point(3,:))/c/dt);
sca=[0 0 0; point]; 
amp=[0; ones(N_sca,1)];

%% Compute STA signals
disp('Computing STA signals');
duration_delay=(2*length(excitation)+1)*dt/2; % delay inserted by Field II 
IR=zeros(cropat,N_elements,N_elements);
t_out=0:dt:((cropat-1)*dt);
wb = waitbar(0, 'Computing IR');
for n=1:N_elements
    waitbar(n/N_elements, wb);

    xdc_apodization(Th,0,[zeros(1,n-1) 1 zeros(1,N_elements-n)]);
    xdc_apodization(Rh,0,ones(1,N_elements));
    [v,t]=calc_scat_multi (Th, Rh, sca, amp);
    t_in=(0:dt:((size(v,1)-1)*dt))+t-duration_delay-2*min_delay; 
    v_aux=interp1(t_in,v,t_out,'linear',0);

    IR(:,:,n)=v_aux;
    
end
close(wb);

% normalising signal
IR=IR./max(abs(IR(:)));

%% Define a reconstruction object
recons=reconstruction();

% define the scan -> only linear scan svailable for the moment 
recons.scan.x_axis=linspace(-20e-3,20e-3,256).';               % x vector [m]
recons.scan.z_axis=linspace(2e-3,42e-3,512).';                 % z vector [m]

% define the transmit & receive beams
%F-number, transmit apodization, steering angle [rad], length of the edge smoothing area [elements], order of the edge smoothing polynomial
recons.transmit_beam=beam(1.2,E.apodization_type.boxcar,0,0,0);
recons.receive_beam=beam(1.2,E.apodization_type.boxcar,0,0,0);

%% Define the sta dataset object
s=sta('Field II, STA, RF format',...    % name of the dataset
      E.signal_format.RF,...            % signal format: RF or IQ
      c,...                             % reference speed of sound (m/s)
      t_out.',...                       % time vector (s)
      IR,...                            % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,2)]);      % probe geometry [x, y, z] (m)

%% reconstruction and show
s.image_reconstruction(recons);
recons.show();

%% demodulate data and create another sta object that avoids the Hilbert artifact
out=tools.demodulate(s);
s2=sta('Field II, STA, IQ format',...     % name of the dataset
      E.signal_format.IQ,...              % signal format: RF or IQ
      c,...                               % reference speed of sound (m/s)
      out.time,...                        % time vector (s)
      out.data,...                        % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,2)],...      % probe geometry [x, y, z] (m)
      out.modulation_frequency);

%% reconstructs and show again
s2.image_reconstruction(recons);
recons.show();
