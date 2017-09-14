function channel_data = create_FI_phased_array_channeldata(h,number_of_frames)

%%%%
%    Save to channeldata for Focused Imaging with phased array imaging
%% Create channel_data object and set some parameters
channel_data = uff.channel_data();
channel_data.sampling_frequency = h.Fs;
channel_data.sound_speed = h.c0;
channel_data.initial_time = 0;
channel_data.probe=create_probe_object(h);

if nargin < 2
    h.number_of_frames = h.Resource.RcvBuffer(1).numFrames;
else
    h.number_of_frames = number_of_frames;
end

%% SEQUENCE GENERATION
N=length(h.TX);                      % number of focused beams
azimuth_axis=h.angles.';
seq=uff.wave();
for n=1:N
    seq(n)=uff.wave();
    seq(n).probe=channel_data.probe;
    
    seq(n).source=uff.point();
    seq(n).source.azimuth=azimuth_axis(n);
    seq(n).source.distance=h.TX(n).focus*h.lambda;
    
    seq(n).apodization = uff.apodization();
    seq(n).apodization.window=uff.window.tukey50;
    seq(n).apodization.f_number=1.7;
    seq(n).apodization.focus=uff.scan('xyz',seq(n).source.xyz);
    
    seq(n).sound_speed=channel_data.sound_speed;
end
channel_data.sequence = seq;

channel_data.pulse = uff.pulse(h.Trans.frequency*10^6);

%% Convert channel data from Verasonics format to USTB format
no_samples = h.Receive(1).endSample;
data = single(zeros(no_samples, channel_data.probe.N, length(seq), h.number_of_frames));
offset_time = calculate_delay_offset(h); % Get offset time
plot_delayed_signal=0;
interpolation_factor = 10;

n=1;
frame_idx = 0;
for n_frame = h.frame_order
    frame_idx = frame_idx + 1;
    for n_tx = 1:length(seq)
        % compute time vector for this line
        t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
        t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
        no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
        
        % compute the offset in time from center of probe to
        % center of transmit wave. We do this by finding the
        % mean between the two center transmit delays for a
        % even numbered probe, and the center transmit delay
        % for a odd elemtn probe
        t0_1 = mean(h.TX(n_tx).Delay(ceil(channel_data.probe.N_elements/2):ceil((channel_data.probe.N_elements+1)/2)))...
            *h.lambda/h.Resource.Parameters.speedOfSound;
        t_out = linspace(t_ini,t_end,no_t);
        t_in=t_out-offset_time-t0_1;
        
        data_in = h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,h.Trans.Connector,h.frame_order(frame_idx));
        data(:,:,n_tx,frame_idx) = time_shift_data(h,data_in,t_in,t_out,interpolation_factor,channel_data);
        n=n+1;
        
        %%
        % to check delay calculation
        % NB! For phased array this is only correct when you
        % are firing at angle=0
        if plot_delayed_signal
            % Point to beamform to (where the scatterer is in the simulation)
            % Need to change to correct scatter setup in the
            % Verasonics script, see FI_phase_array_p4.m for
            % example. This seems to be correct, but the delays
            % are slighty off for transmit angles > 0 but not
            % for angles < 0. Strange. Is there somthing wrong
            % with the Verasonics simulation?? :)
            
            [z_all,x_all] = pol2cart(h.angles,ones(1,length(h.angles))*40e-3);
            x = x_all(n_tx);
            y = 0;
            z = z_all(n_tx);
            
            TF=(-1).^(z<channel_data.sequence(n_tx).source.z).*sqrt((channel_data.sequence(n_tx).source.x-x).^2+(channel_data.sequence(n_tx).source.y-y).^2+(channel_data.sequence(n_tx).source.z-z).^2);
            % add distance from source to origin
            TF=TF+sign(cos(channel_data.sequence(n_tx).source.azimuth)).*channel_data.sequence(n_tx).source.distance;
            % receive delay
            RF=sqrt((channel_data.probe.x-x).^2+(channel_data.probe.y-y).^2+(channel_data.probe.z-z).^2);
            % total delay
            delay=(RF+TF)/channel_data.sound_speed;
            %%
            figure(101); hold off;
            pcolor(1:channel_data.probe.N_elements,time,real(data(:,:,n_tx,n_frame))); shading flat; colormap gray; colorbar; hold on;
            plot(1:channel_data.probe.N_elements,delay,'r');
            title(n_tx);
            ylim([0.9*min(delay) 1.1*max(delay)]);
            pause();
        end
    end
end

channel_data.data = data;

%%

end