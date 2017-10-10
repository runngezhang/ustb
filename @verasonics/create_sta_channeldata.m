function channel_data = create_sta_channeldata(h)
%%%%
%    Save to channeldata for Synthetic Transmit Aperture imaging

%% Create channel_data object and set some parameters
channel_data = uff.channel_data();
channel_data.sampling_frequency = h.Fs;
channel_data.sound_speed = h.c0;
channel_data.initial_time = 0;
channel_data.probe=create_probe_object(h);


%% SEQUENCE GENERATION
N=length(h.TX);                      % number of waves
for n=1:N
    seq(n)=uff.wave();
    seq(n).probe=channel_data.probe;
    seq(n).source.xyz=[channel_data.probe.x(n) channel_data.probe.y(n) channel_data.probe.z(n)];
    
    seq(n).apodization = uff.apodization();
    seq(n).apodization.window=uff.window.sta;
    seq(n).apodization.origo=uff.point('xyz', [0 0 -Inf]);
    
    seq(n).sound_speed=channel_data.sound_speed;
end
channel_data.sequence = seq;

% Add center frequency to channel_data
channel_data.pulse = uff.pulse();
channel_data.pulse.center_frequency = h.Trans.frequency*10^6;

%% Convert channel data from Verasonics format to USTB format
no_samples = h.Receive(1).endSample;
data = single(zeros(no_samples, h.Resource.Parameters.numRcvChannels, length(seq), h.Resource.RcvBuffer(1).numFrames));

offset_time = calculate_delay_offset(h); % Get offset time
%time = [0:(1/h.Fs):((no_samples-1)/h.Fs)]';
plot_delayed_signal=0;
n=1;
interpolation_factor = 10;
for n_frame = h.frame_order
    for n_tx = 1:length(seq)
        % compute time vector for this line
        t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
        t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
        no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
        
        % compute the offset in time from center of probe to
        % element to get correct t_0
        element_offset = channel_data.probe.r(n_tx)/channel_data.sound_speed;
        
        t_out = linspace(t_ini,t_end,no_t);
        t_in=t_out-offset_time+element_offset;
        
        % read data
        data_in = double(h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,:,n_frame));
        data(:,:,n_tx,n_frame) = time_shift_data(h,data_in,t_in,t_out,interpolation_factor,channel_data);
        %data(:,:,n_tx,n_frame)=interp1(t_in,,time,'linear',0);
        n=n+1;
        
        
        % to check delay calculation
        if plot_delayed_signal
            % Point to beamform to (where the scatterer is in the simulation)
            x = 0;
            y = 0;
            z = 20e-3;
            
            TF=(-1).^(z<channel_data.sequence(n_tx).source.z).*sqrt((channel_data.sequence(n_tx).source.x-x).^2+(channel_data.sequence(n_tx).source.y-y).^2+(channel_data.sequence(n_tx).source.z-z).^2);
            % add distance from source to origin
            TF=TF+sign(cos(channel_data.sequence(n_tx).source.azimuth)).*channel_data.sequence(n_tx).source.distance;
            % receive delay
            RF=sqrt((channel_data.probe.x-x).^2+(channel_data.probe.y-y).^2+(channel_data.probe.z-z).^2);
            % total delay
            delay=(RF+TF)/channel_data.sound_speed;
            
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

end
