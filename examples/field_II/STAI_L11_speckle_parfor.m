%% Computation of a STAI dataset with Field II and beamforming with USTB
%
% This example shows how to load the data from a Field II simulation into 
% USTB objects, and then beamformt it with the USTB routines. 
% This example uses the L11-4v 128 element Verasonics Transducer
% The Field II simulation program (field-ii.dk) should be in MATLAB's path.
%
% date:     11.03.2015
% updated:  09.05.2017
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%           Ole Marius Hoel Rindal <olemarius@olemarius.net>

clear all;
close all;

%% basic constants
c0=1540;     % Speed of sound [m/s]
fs=100e6;    % Sampling frequency [Hz]
dt=1/fs;     % Sampling step [s] 

%% field II initialisation
field_init;
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% Transducer definition L11-4v, 128-element linear array transducer
probe = uff.linear_array();
f0                      = 5.1333e+06;      % Transducer center frequency [Hz]
lambda                  = c0/f0;           % Wavelength [m]
probe.element_height    = 5e-3;            % Height of element [m]
probe.pitch             = 0.300e-3;        % probe.pitch [m]
kerf                    = 0.03e-03;        % gap between elements [m]
probe.element_width     = probe.pitch-kerf;% Width of element [m]
lens_el                 = 20e-3;           % position of the elevation focus
probe.N                 = 128;             % Number of elements
pulse_duration          = 2.5;             % pulse duration [cycles]

%% pulse definition
pulse = uff.pulse(f0);
pulse.fractional_bandwidth = 0.65;        % probe bandwidth [1]
t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);
impulse_response = gauspuls(t0, f0, pulse.fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
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
noSubAz=round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl=round(probe.element_height/(lambda/8));       % number of subelements in the elevation direction
Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 

% setting excitation, impulse response and baffle
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);

%% phantom of speckle
number_of_scatterers = 100;
xxp_speckle=random('unif',-5e-3,5e-3,number_of_scatterers,1);
zzp_speckle=random('unif',15e-3,20e-3,number_of_scatterers,1);
sca = [xxp_speckle zeros(length(xxp_speckle),1) zzp_speckle];  % list with the scatterers coordinates [m]
amp=randn(length(sca));                   % list with the scatterers amplitudes
cropat=round(1.1*2*sqrt((max(sca(:,1))-min(probe.x))^2+max(sca(:,3))^2)/c0/dt);   % maximum time sample, samples after this will be dumped
%% output data
t_out=0:dt:((cropat-1)*dt);                 % output time vector
STA=zeros(cropat,probe.N,probe.N);    % impulse response channel data
%% Compute STA signals
disp('Field II: Computing STA dataset');
disp('No waitbar possible for parfor, so just be patient :)');
parfor n=1:probe.N
    %Since we are using parfor, we have to initate Field II and the arrays
    %for every worker as well.
    field_init(0);
    Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
    Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);
    xdc_excitation (Th, excitation);
    xdc_impulse (Th, impulse_response);
    xdc_baffle(Th, 0);
    xdc_center_focus(Th,[0 0 0]);
    xdc_impulse (Rh, impulse_response);
    xdc_baffle(Rh, 0);
    xdc_center_focus(Rh,[0 0 0]);
    
    % transmit aperture
    xdc_apodization(Th, 0, [zeros(1,n-1) 1 zeros(1,probe.N-n)]);
    xdc_focus_times(Th, 0, zeros(1,probe.N));
    
    % receive aperture    
    xdc_apodization(Rh, 0, ones(1,probe.N));
    xdc_focus_times(Rh, 0, zeros(1,probe.N));
    
    % do calculation
    [v,t]=calc_scat_multi(Th, Rh, sca, amp);
    
    % lag compensation
    t_in=(0:dt:((size(v,1)-1)*dt))+t-lag*dt + probe.r(n)/c0;
    v_aux=interp1(t_in,v,t_out,'linear',0);

    % build the dataset
    STA(:,:,n)=v_aux;
    
    %% SEQUENCE GENERATION
    seq(n)=uff.wave();
    seq(n).probe=probe;
    seq(n).source.xyz=[probe.x(n) probe.y(n) probe.z(n)];
    seq(n).sound_speed=c0;
end

%% CHANNEL DATA
channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.data = STA*10^29;

%% SCAN
scan=uff.linear_scan(linspace(-5e-3,5e-3,256).', linspace(15e-3,20e-3,256).');

%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=scan;
bmf.receive_apodization.window=uff.window.boxcar;
bmf.receive_apodization.f_number=1.7;
bmf.receive_apodization.apex.distance=Inf;
bmf.transmit_apodization.window=uff.window.boxcar;
bmf.transmit_apodization.f_number=1.7;
bmf.transmit_apodization.apex.distance=Inf;

% Delay and sum on receive, then coherent compounding
b_data=bmf.go({process.das_matlab() process.coherent_compounding()});
% Display image
b_data.plot()

%%
envelope = abs(b_data.data);
envelope = envelope./max(envelope(:));
m = mean(envelope(:));
s = std(envelope(:));

snr_calculated_das = m/s
snr_theoretical = (pi/(4-pi))^(1/2)
b = s/(sqrt((4-pi)/2)); %Scale parameter

% Estimate PDF
x_axis = linspace(0,1,200);
[n,xout] = hist(envelope(:),x_axis);
delta_x = xout(2)-xout(1);
n = n/sum(n)/delta_x;

% Theoretical Raileigh PDF 
theoretical_pdf = (x_axis./b^2).*exp(-x_axis.^2/(2.*b^2));

% Plot
color=[0.25 1 0.75]
figure(1);clf;  
plot(xout,n,'LineWidth',2,'Color','r','DisplayName','Estimated PDF');hold on;
plot(x_axis,theoretical_pdf,'--','Color',color,'LineWidth',2,'DisplayName','Rayleigh Theoretical PDF');
title('PDF of envelope');
xlabel('Normalized amplitude');
ylabel('Probability')
legend('show');

%% Save UFF dataset
uff_file=uff('FieldII_speckle_simulation.uff');
uff_file.write(channel_data,'channel_data');
uff_file.write(b_data,'b_data');