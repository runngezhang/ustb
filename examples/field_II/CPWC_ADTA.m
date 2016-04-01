%% Implementation of Angular Dependent Transmit Apodization for a CPWC 
%
%
% created:  26.03.2016
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

clear all;
close all;

%% basic constants
c0=1540;     % Speed of sound [m/s]
fs=100e6;   % Sampling frequency [Hz]
dt=1/fs;    % Sampling step [s] 

%% field II initialisation
field_init;
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% transducer definition L11 Ultrasonix
f0=5.2e6;           % Transducer center frequency [Hz]
bw=0.67;            % probe bandwidth [1]
lambda=c0/f0;       % Wavelength [m]
height=5e-3;        % Height of element [m]
pitch=300e-6;       % pitch [m]
kerf=30e-6;      % gap between elements [m]
width=pitch-kerf;   % Width of element [m]
lens_el=20e-3;      % position of the elevation focus
N_elements=128;     % Number of elements
cycles=1.5;

%% pulse definition
t0 = (-1/bw/f0): dt : (1/bw/f0);
impulse_response = gauspuls(t0, f0, bw);
te = (-cycles/2/f0): dt : (cycles/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2;   

% show the pulse to check that the lag estimation is on place (and that the pulse is symmetric)
figure;
plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
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
Na=3;                                       % number of plane waves 
angles=linspace(-alpha_max,alpha_max,Na);   % vector of angles [rad]
F=1;                                        % number of frames
%angles=repmat(angles,[1 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADA (Angle Dependent Apodization)
P=cycles/f0;     % pulse duration


% writing sequence angles
apo=zeros(length(angles),N_elements);
for n = 1:length(angles);
    switch round((1+n)/3)
        case 1
            xc=-10e-3; zc=20e-3;     % sweet point (center of the image)
        case 2
            xc=0; zc=20e-3;     % sweet point (center of the image)
        case 3
            xc=10e-3; zc=20e-3;     % sweet point (center of the image)
    end
    % lateral distance
    xd=x0-(xc-zc*tan(angles(n)));
    % theoretical limit
    Dxp=(sqrt(P^2*c0^2+4*zc*P*c0*cos(angles(n)))-P*c0*sin(angles(n)))/2/cos(angles(n))^2;
    Dxn=(sqrt(P^2*c0^2+4*zc*P*c0*cos(angles(n)))+P*c0*sin(angles(n)))/2/cos(angles(n))^2;
    % build apodization            
    Apod=zeros(size(xd));
    Apod(xd>=0)=abs(xd(xd>=0))<Dxp;
    Apod(xd<0)=abs(xd(xd<0))<Dxn;
    % smooth it a bit -> not really needed
    fApod=filtfilt(ones(1,5)/5,1,Apod);
    apo(n,:)=fApod;

    figure(101);
    plot(Apod); hold on; grid on;
    plot(fApod,'r'); 
end
%% end of ADA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% phantom
x_point=0;%-15e-3:1e-4:15e-3;                                      % x-coordinate of the scatterers [m]
z_point=15e-3;%[13e-3:1e-4:14e-3 22e-3:1e-4:23e-3];%5e-3:5e-3:40e-3;   % z-coordinate of the scatterers [m]
[xxp zzp]=meshgrid(x_point,z_point);
N_sca=length(zzp(:));                                       % total number of scatterers
sca=[xxp(:) zeros(N_sca,1) zzp(:)];                         % list with the scatterers coordinates [m]
amp=ones(N_sca,1);                                          % list with the scatterers amplitudes
cropat=round(2*sqrt((max(x_point)-min(x_point))^2+max(z_point)^2)/c0/dt);   % maximum time sample, samples after this will be dumped
amp=random('unif',0,1,size(sca,1),1);

%% output data
t_out=0:dt:((cropat-1)*dt);         % output time vector
CPW=zeros(cropat,N_elements,length(angles),F);  % impulse response channel data

%% Compute CPW signals
time_index=0;
disp('Field II: Computing CPW dataset');
wb = waitbar(0, 'Field II: Computing CPW dataset');
for f=1:F
    waitbar(0, wb, sprintf('Field II: Computing CPW dataset, frame %d',f));
    for n=1:length(angles)
        waitbar(n/length(angles), wb);
        
        % transmit aperture
        %xdc_apodization(Th,0,apo(n,:));
        xdc_apodization(Th,0,ones(1,N_elements));
        xdc_times_focus(Th,0,x0(:).'*sin(angles(n))/c0);

        % receive aperture
        xdc_apodization(Rh, 0, ones(1,N_elements));
        xdc_focus_times(Rh, 0, zeros(1,N_elements));
        
        % do calculation
        [v,t]=calc_scat_multi(Th, Rh, sca, amp);

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
      c0,...                             % reference speed of sound (m/s)
      angles.',...                       % angle vector [rad]
      t_out.',...                       % time vector (s)
      CPW,...                            % matrix with the data [samples, channels, firings, frames]
      [x0.' zeros(N_elements,2)]);      % probe geometry [x, y, z] (m)

% convert to IQ data
cpw_dataset.demodulate(true,[],[],[],E.demodulation_algorithm.oyvind);

%data_1=cpw_dataset.data(:,:,1:3);
%data_2=cpw_dataset.data(:,:,4:6);
%data_3=cpw_dataset.data(:,:,7:9);

%% define a reconstruction

% define a scan
scan=linear_scan();
scan.x_axis=linspace(-15e-3,15e-3,1024).';               % x vector [m]
scan.z_axis=linspace(7e-3,27e-3,512).';                 % z vector [m]

% define a synthetic orientation
F_number=1.7;
orien=orientation();
orien.transmit_beam=beam(F_number,E.apodization_type.none);  
orien.receive_beam=beam(F_number,E.apodization_type.tukey25);

% define a reconstruction 
cpw_image=reconstruction();
cpw_image.scan=scan;
cpw_image.orientation=orien;

%% beamform and show
cpw_dataset.angle=cpw_dataset.angle(1:3)

cpw_dataset.data=data_1;
cpw_dataset.image_reconstruction(cpw_image);
image_1=cpw_image.show();

cpw_dataset.data=data_2;
cpw_dataset.image_reconstruction(cpw_image);
image_2=cpw_image.show();

cpw_dataset.data=data_3;
cpw_dataset.image_reconstruction(cpw_image);
image_3=cpw_image.show();


mask_1=abs(scan.x_axis+10e-3)<6e-3; mask_1=filtfilt(ones(1,40)/40,1,double(mask_1)); mask_1=mask_1*60-60;
mask_2=(abs(scan.x_axis+0e-3)<6e-3); mask_2=filtfilt(ones(1,40)/40,1,double(mask_2)); mask_2=mask_2*60-60;
mask_3=(abs(scan.x_axis-10e-3)<6e-3); mask_3=filtfilt(ones(1,40)/40,1,double(mask_3)); mask_3=mask_3*60-60;
figure;
plot(mask_1); hold on;
plot(mask_2);
plot(mask_3);
final_image=log10(bsxfun(@times,10.^image_1,10.^mask_1.')+bsxfun(@times,10.^image_2,10.^mask_2.')+bsxfun(@times,10.^image_3,10.^mask_3.'));
final_image=image_1;

figure;
imagesc(scan.x_axis*1e3,scan.z_axis*1e3,final_image);
colormap gray; caxis([-50 0]); 
colorbar;
axis equal tight;  
xlabel('x[mm]');
ylabel('z[mm]');
set(gca,'fontsize', 18);
title('FLAT');
drawnow;

