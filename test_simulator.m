%% PHANTOM
pha=phantom();
pha.sound_speed=1540;               % speed of sound [m/s]
pha.points=[0,  0,  0, 0;
                0,  0, 20e-3, 1;    % point scatterer position [m]
            -5e-3,  0, 20e-3, 1; 
             5e-3,  0, 20e-3, 1; 
            -10e-3, 0, 20e-3, 1;
             10e-3, 0, 20e-3, 1;
                 0, 0,  5e-3, 1;
                 0, 0, 10e-3, 1;
                 0, 0, 15e-3, 1;
                 0, 0, 25e-3, 1;
                 0, 0, 30e-3, 1;
                 0, 0, 35e-3, 1];

% checking phantom
pha.plot();             
             
%% PROBE
prob=probe();
N_elements=128;                              % number of elements in the array
pitch=300e-6;                                % probe pitch
a=270e-6;                                    % element width [m]
b=5e-3;                                      % element height [m]
x0=((1:N_elements)*pitch).'; x0=x0-mean(x0); % element position in the x_axis (m)
y0=zeros(N_elements,1);                      % element position in the y_axis [m]
z0=zeros(N_elements,1);                      % element position in the z_axis [m]
theta=zeros(N_elements,1);                   % element orientation in the azimuth direction [rad]
phi=zeros(N_elements,1);                     % element orientation in the elevation direction [rad]
prob.geometry=[x0 y0 z0 theta phi a*ones(N_elements,1) b*ones(N_elements,1)]; % probe geometry
  
% checking probe
prob.plot();

%% Pulse
%
% The pulse class implements a Gaussian-modulated RF pulse of known center
% frequency and fractional bandsidth.

tx_pulse=pulse();
tx_pulse.center_frequency=5.2e6;                           % transducer frequency [MHz]
tx_pulse.fractional_bandwidth=0.6;                         % fractional bandwidth [unitless]

% checking pulse
tx_pulse.plot();

%% The beam sequence
sou=source();
source.xyz=


N_active_elements=32;                                   % number of active elements
number_beams = N_elements - (N_active_elements-1);      % number of beams
focus_length=20e-3;                                     % focus length

sequence=repmat(beam(),1,number_beams);
for n_beam=1:number_beams 
    
    % select active channels
    sequence(n_beam).apodization=[zeros(1,n_beam-1) ones(1,N_active_elements) zeros(1,N_elements-(n_beam+N_active_elements)+1,1)].';
  
    % compute transmit delays
    beam_center_x=prob.x.'*sequence(n_beam).apodization/sum(sequence(n_beam).apodization);
    sequence(n_beam).delay=(focus_length-sqrt(sum((prob.geometry(:,1:3)-ones(N_elements,1)*[beam_center_x 0 focus_length]).^2,2)))/pha.sound_speed;
    
    % show
    sequence(n_beam).plot();
end


%% The model
%
% Here comes the model.

my_model=model();

% setting input data 
my_model.phantom=pha;               % phantom
my_model.pulse=tx_pulse;                    % transmitted pulse
my_model.probe=prob;                % probe
my_model.sequence=sequence;                 % beam sequence
my_model.sampling_frequency=41.6e6;        % sampling frequency [Hz]

% we launch the simulation
my_dataset=my_model.simulate();

% check how does it look
my_dataset.plot(48);

%% Beamforming the dataset
%
% This is not part of this assignment, but it would be silly not to check
% the image

% DRF beamforming
sr_image=zeros(my_dataset.N_samples,my_dataset.N_beams);
x_axis=zeros(my_dataset.N_beams,1);
for n_beam=1:my_dataset.N_beams
    x_axis(n_beam)=sum(my_dataset.sequence(n_beam).apodization.*prob.x)/sum(my_dataset.sequence(n_beam).apodization);
    % computing the receive signal
    for n_rx=1:my_dataset.N_elements
        if(sequence(n_beam).apodization(n_rx)>0)
            focus_vector=[x_axis(n_beam)*ones(my_dataset.N_samples,1) zeros(my_dataset.N_samples,1) my_dataset.phantom.sound_speed*my_dataset.time.'/2];
            rx_focusing_delay=(focus_vector(:,3)-sqrt(sum((ones(size(my_dataset.time,2),1)*prob.geometry(n_rx,1:3)-focus_vector).^2,2)))/pha.sound_speed;  % dynamic focussing delay
            sr_image(:,n_beam)=sr_image(:,n_beam)+interp1(my_dataset.time,my_dataset.data(:,n_rx,n_beam),my_dataset.time-rx_focusing_delay.','linear',0).';           % signal
        end
    end
end

% convert to intensity values
envelope_drf=abs(hilbert(sr_image));
envelope_drf_dB=20*log10(envelope_drf./max(envelope_drf(:)));

figure3 = figure('Color',[1 1 1]); 
imagesc(x_axis*1e3,my_dataset.time*pha.sound_speed/2*1e3,envelope_drf_dB); axis tight equal; 
box('on'); 
xlabel('x [mm]');
ylabel('z [mm]')
set(figure3,'InvertHardcopy','off');
caxis([-60 0]); colorbar; colormap gray;
set(gca,'color','black')
ylim([0 40]);

