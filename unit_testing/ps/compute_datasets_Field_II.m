clear all;
close all;

% field initialisation
field_init;
set_field('att', 0);
c=1540; set_field('c',c);                % Speed of sound [m/s]
fs=100e6; set_field('fs',fs); dt=1/fs;    % Sampling frequency [Hz]
set_field('use_rectangles',1);

% Define transducer
f0=4.5e6;           % Transducer center frequency [Hz]
lambda=c/f0;        % Wavelength [m]
height=4e-3;        % Height of element [m]
pitch=0.33e-3;      % pitch [m]
kerf=pitch./10;     % Distance between transducer elements [m]
width=pitch-kerf;   % Width of element [m]
N_elements=128;     % Number of elements

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

% setting a minimal delay 
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
point=[0e-3 0e-3 40e-3];
cropat=round(2*2*point(3)/c/dt);
N_sca=1;
sca=[0 0 0; point; 0 0 2.5*point(3)]; % single scatterer
amp=[0; 1; 0];

% Compute IR signals
disp('Computing STA signals');
duration_delay=(2*length(excitation)+1)*dt/2; % delay inserted by Field II 
IR=[];
t_out=0:dt:((cropat-1)*dt);
%times=zeros(1,N_elements);
wb = waitbar(0, 'Computing IR');
for n=1:N_elements
    waitbar(n/N_elements, wb);

    xdc_apodization(Th,0,[zeros(1,n-1) 1 zeros(1,N_elements-n)]);
    xdc_apodization(Rh,0,ones(1,N_elements));
    [v,t]=calc_scat_multi (Th, Rh, sca, amp);
    %times(n)=t;
    t_in=(0:dt:((size(v,1)-1)*dt))+t-duration_delay-2*min_delay; 
    v_aux=interp1(t_in,v,t_out,'linear',0);

    IR(:,:,n)=v_aux;
    
end
close(wb);

% normalising signal
IR=IR./max(abs(IR(:)));

%% STA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving STA-IR dataset
disp('Saving STA signals');
s=struct(   'name','Point scatterer (sta) RF',...   % name of the dataset
            'format',E.signal_format.RF,...         % signal format: RF or IQ
            'c0',c,...                              % reference speed of sound (m/s)
            'time',t_out.',...                        % time vector (s)
            'data',IR,...                           % matrix with the data [samples, channels, firings, frames]
            'geom',[x0.' zeros(N_elements,2)]);     % probe geometry [x, y, z] (m)
save('./struct/ps_sta_rf.mat','s');

% saving STA-IQ dataset
s=tools.demodulate(s);
s.name='Point scatterer (sta) IQ';
save('./struct/ps_sta_iq.mat','s');

%% CPW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing CPW signals');

% chosing angle squence 
FN_min=1.75;
alpha_max=atan(1/2/FN_min);
lambda=c/f0;
Na=71; 
alpha=linspace(-alpha_max,alpha_max,Na);

% computing CPW dataset from STA dataset
PW=zeros(size(IR,1),N_elements,Na);
wb=waitbar(0,'CPW from STA');
for na=1:Na
    waitbar(na/length(alpha));
    for tx=1:N_elements
        PW(:,:,na)= PW(:,:,na)+interp1(t_out,IR(:,:,tx),t_out-x0(tx)*sin(alpha(na))/c,'linear',0);
    end
end
close(wb);

% saving CPW-RF dataset
disp('Saving CPW signals');
s=struct(   'name','Point scatterer (sta) RF',...   % name of the dataset
            'format',E.signal_format.RF,...         % signal format: RF or IQ
            'c0',c,...                              % reference speed of sound (m/s)
            'angle',alpha.',...                     % vector of angles (rad)
            'time',t_out.',...                        % time vector (s)
            'data',PW,...                           % matrix with the data [samples, channels, firings, frames]
            'geom',[x0.' zeros(N_elements,2)]);     % probe geometry [x, y, z] (m)
   
save('./struct/ps_cpw_rf.mat','s');   

% saving CPW-IQ dataset
s=tools.demodulate(s);
s.name='Point scatterer (cpw) IQ';
save('./struct/ps_cpw_iq.mat','s');

%% VS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% chosing position of virtual sources
Na=70; 
vx=linspace(min(x0),max(x0),Na);
vz=-10e-3.*ones(1,Na);

% computing VS dataset from STA
VS=zeros(size(IR,1),N_elements,Na);
wb=waitbar(0,'VS from STA');
for na=1:Na
    waitbar(na/Na);
    for tx=1:N_elements
        delays=sqrt((x0(tx)-vx(na)).^2+vz(na).^2)/c;
        VS(:,:,na)= VS(:,:,na)+interp1(t_out,IR(:,:,tx),t_out-delays,'linear',0);
    end
end
close(wb);

% saving VS-RF dataset
s=struct(   'name','Point scatterer (vs) RF',...    % name of the dataset
            'format',E.signal_format.RF,...         % signal format: RF or IQ
            'c0',c,...                              % reference speed of sound (m/s)
            'source',[vx.' zeros(Na,1) vz.'],...    % position of the virtual sources [x, y, z] (m)
            'time',t_out.',...                      % time vector (s)
            'data',VS,...                           % matrix with the data [samples, channels, firings, frames]
            'geom',[x0.' zeros(N_elements,2)]);     % probe geometry [x, y, z] (m)
save('./struct/ps_vs_rf.mat','s');

% saving VS-IQ dataset
s=tools.demodulate(s);
s.name='Point scatterer (vs) IQ';
save('./struct/ps_vs_iq.mat','s');
