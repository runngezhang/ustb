function channel_data = create_cpw_channeldata(h)
%% Create channel_data object and set some parameters
channel_data = uff.channel_data();
channel_data.sampling_frequency = h.Fs;
channel_data.sound_speed = h.c0;
channel_data.initial_time = 0;
channel_data.probe=create_probe_object(h);

%% SEQUENCE GENERATION
N=size(h.TX,2);             % number of plane waves
for n=1:N
    seq(n)=uff.wave();
    seq(n).probe=channel_data.probe;
    seq(n).source.azimuth=h.angles(n);
    seq(n).source.distance=Inf;
    seq(n).sound_speed=channel_data.sound_speed;
end
channel_data.sequence = seq;


%% Save Pulse
channel_data.pulse = uff.pulse(double(h.Trans.frequency*10^6));

%% Convert channel data from Verasonics format to USTB format
no_samples = h.Receive(1).endSample;
data = zeros(no_samples, h.Resource.Parameters.numRcvChannels, length(seq), h.Resource.RcvBuffer(1).numFrames);

offset_time = calculate_delay_offset(h); % Get offset time
n=1;
time = [0:(1/h.Fs):((no_samples-1)/h.Fs)]';
plot_delayed_signal=0;
for n_frame = h.frame_order
    for n_tx = 1:length(channel_data.sequence)
        %% compute time vector for this line
        t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
        t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
        no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
        
        % Find t_0, when the plane wave "crosses" the center of
        % the probe
        if 1  %Calculate geometrically
            D = abs(h.Trans.ElementPos(1,1)-h.Trans.ElementPos(end,1))*1e-3;
            q = abs((D/2)*sin(channel_data.sequence(n_tx).source.azimuth));
            t0_1 = q/(channel_data.sound_speed);
        else  %Calculate using Verasonics transmit delay, this will not work for the multiplexer probe NBNB!
            t0_1 = mean(h.TX(n_tx).Delay)*h.lambda/h.Resource.Parameters.speedOfSound;
            figure(100);hold all;
            plot(h.TX(n_tx).Delay)
            plot((h.TX(n_tx).Delay(end/2))*ones(1,128),'r')
            plot(mean(h.TX(n_tx).Delay)*ones(1,128),'b')
        end
        
        t_in=linspace(t_ini,t_end,no_t)-offset_time-t0_1;
        
        if isfield(h.Trans,'HVMux') % If Transducer has MUX, for example the L12-4v, we need to re arrange channels
            validChannels = h.Trans.HVMux.Aperture(:,h.Receive(n_tx).aperture)';
            validChannels = validChannels(validChannels>0);
        else
            validChannels = [1:128];
        end
        %% read data
        data(:,:,n_tx,n_frame)=interp1(t_in,double(h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,validChannels,n_frame)),time,'linear',0);
        n=n+1;
        %%
        % to check delay calculation
        if plot_delayed_signal
            %delay= 20e-3*cos(angles(n_tx))/h.c0+delay_x0;
            %%
            z = 20e-3;
            x = 0;
            y = 0;
            TF = z*cos(channel_data.sequence(n_tx).source.azimuth)*cos(channel_data.sequence(n_tx).source.elevation)+x*sin(channel_data.sequence(n_tx).source.azimuth)*cos(channel_data.sequence(n_tx).source.elevation)
            % receive delay
            RF=sqrt((channel_data.probe.x-x).^2+(channel_data.probe.y-y).^2+(channel_data.probe.z-z).^2);
            % total delay
            delay=(RF+TF)/channel_data.sound_speed;
            
            figure(101); hold off;
            pcolor(1:length(channel_data.probe.x),time,abs(data(:,:,n_tx,n_frame))); shading flat; colormap gray; colorbar; hold on;
            plot(1:length(channel_data.probe.x),delay,'r');
            title(n_tx);
            ylim([0.95*min(delay) 1.05*max(delay)]);
            pause();
        end
    end
end

channel_data.data = data;

end