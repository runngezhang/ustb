classdef verasonics < handle
    %Verasonics  Class for reading Verasonics research scanner data to USTB
    
    %   authors:    Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %               Alfonso Rodriques-Morales <alfonso.r.molares@ntnu.no>
    %
    %   $Date: 2017/03/16$
    
    properties (SetAccess = public)
        % Verasonics objects and structs
        Trans                  % Verasonics Transducer object
        TW                     % Verasonics Transmit Waveform object
        RcvData                % Verasonics Receive Data buffer, the channeldata
        Resource               % Verasonics Resource object, a few global parameters
        Receive                % Verasonics Receive object, defining receive parameters
        TX                     % Verasonics TX object, defining transmit parameters
        
        % Some helpful parameters
        Fs                     % The sampling frequency in Hz
        f0                     % The center frequency in Hz
        c0                     % The speed of sound in m/s
        lambda                 % The wavelength in m
        
        % For CPW
        angles                 % The transmit angles of the planewaves in radians, or for pha the transmit beam angles
    end
    
    %% Constructor
    methods (Access = public)
        function h = verasonics()
            %empty constructor
        end
    end
    
    %% Set methods
    
    methods
        function set.Trans(h,Trans)
            assert(strcmp(Trans.units,'mm'),'Please use mm as units in Verasonics.');
            h.Trans = Trans;
            h.f0 = Trans.frequency*10^6;
        end
        
        function set.Receive(h,Receive)
            assert(isempty(h.Trans)==0,'Please set the Trans variable first.');
            h.Receive = Receive;
            h.Fs = h.f0*Receive(1).samplesPerWave;
            if isfield(Receive,'aperture') == 0 % Then this is a no-mux probe and we set this to one
                h.Receive(1).aperture = 1;
            end
        end
        
        function set.Resource(h,Resource)
            assert(isempty(h.Trans)==0,'Please set the Trans variable first.');
            h.Resource = Resource;
            h.c0 = Resource.Parameters.speedOfSound;
            h.lambda = h.c0/h.f0;
        end
        function set.TW(h,TW)
            assert(isempty(h.Trans)==0,'Please set the Trans variable first.');
            h.TW=TW;
        end
    end
    
    %% Create datasets
    methods (Access = public)
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
            for n_frame = 1:h.Resource.RcvBuffer(1).numFrames
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
        
        
        %%%%
        %    Save to channeldata for Synthetic Transmit Aperture imaging
        function channel_data = create_sta_channeldata(h)
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
                seq(n).apodization.apex=seq(n).source;
                seq(n).sound_speed=channel_data.sound_speed;
            end
            channel_data.sequence = seq;
            
            % Add center frequency to channel_data
            channel_data.pulse = uff.pulse(h.Trans.frequency*10^6);
            
            %% Convert channel data from Verasonics format to USTB format
            no_samples = h.Receive(1).endSample;
            data = zeros(no_samples, h.Resource.Parameters.numRcvChannels, length(seq), h.Resource.RcvBuffer(1).numFrames);
            
            offset_time = calculate_delay_offset(h); % Get offset time
            time = [0:(1/h.Fs):((no_samples-1)/h.Fs)]';
            plot_delayed_signal=0;
            n=1;
            for n_frame = 1:h.Resource.RcvBuffer(1).numFrames
                for n_tx = 1:length(seq)
                    % compute time vector for this line
                    t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
                    t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
                    no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
                    
                    % compute the offset in time from center of probe to
                    % element to get correct t_0
                    element_offset = channel_data.probe.r(n_tx)/channel_data.sound_speed;
                    
                    t_in=linspace(t_ini,t_end,no_t)-offset_time+element_offset;
                    
                    % read data
                    data(:,:,n_tx,n_frame)=interp1(t_in,double(h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,:,n_frame)),time,'linear',0);
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
        
        %%%%
        %    Save to channeldata for Focused Imaging with phased array imaging
        function channel_data = create_FI_phased_array_channeldata(h,number_of_frames)
            %% Create channel_data object and set some parameters
            channel_data = uff.channel_data();
            channel_data.sampling_frequency = h.Fs;
            channel_data.sound_speed = h.c0;
            channel_data.initial_time = 0;
            channel_data.probe=create_probe_object(h);
            
            if nargin < 2
                number_of_frames = h.Resource.RcvBuffer(1).numFrames;
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
                seq(n).apodization.scan.xyz=seq(n).source.xyz;
                
                seq(n).sound_speed=channel_data.sound_speed;
            end
            channel_data.sequence = seq;
            
            channel_data.pulse = uff.pulse(h.Trans.frequency*10^6);
            
            %% Convert channel data from Verasonics format to USTB format
            no_samples = h.Receive(1).endSample;
            data = (zeros(no_samples, channel_data.probe.N, length(seq), number_of_frames));
            offset_time = calculate_delay_offset(h); % Get offset time
            plot_delayed_signal=0;
            interpolation_factor = 10;
        
            n=1;
            for n_frame = 1:number_of_frames
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

                    data_in = h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,h.Trans.Connector,n_frame);
                    data(:,:,n_tx,n_frame) = time_shift_data(h,data_in,t_in,t_out,interpolation_factor,channel_data);
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
        
        %%%%
        %    Save to channeldata for Focused Imaging with phased array imaging
        function channel_data = create_FI_linear_array_channeldata(h)
            %% Create channel_data object and set some parameters
            channel_data = uff.channel_data();
            channel_data.sampling_frequency = h.Fs;
            channel_data.sound_speed = h.c0;
            channel_data.initial_time = 0;
            channel_data.probe=create_probe_object(h);
            
            if strcmp(h.TW.type,'parametric') % read pulse fr. from TW
                channel_data.pulse=uff.pulse(h.TW.Parameters(1)*1e6,0,0);
            else % read pulse fr. from transducer
                channel_data.pulse=uff.pulse(h.f0,0,0);
            end
            
            
            %% SEQUENCE GENERATION
            N=length(h.TX);                 % number of focused beams
            seq=uff.wave();
            for n=1:N
                seq(n)=uff.wave();
                seq(n).probe=channel_data.probe;
                seq(n).source.xyz=[h.TX(n).Origin(1)*h.lambda 0 h.TX(n).focus*h.lambda];
                
                seq(n).apodization = uff.apodization();
                seq(n).apodization.window=uff.window.tukey50;
                seq(n).apodization.f_number=1.7;
                seq(n).apodization.apex.distance=Inf;
                seq(n).apodization.scan.xyz=seq(n).source.xyz;
                
                seq(n).sound_speed=channel_data.sound_speed;
                
                % show source
            end
            channel_data.sequence = seq;
            %% Convert channel data from Verasonics format to USTB format
            no_samples = h.Receive(1).endSample;
            data = zeros(no_samples, channel_data.probe.N, length(seq), h.Resource.RcvBuffer(1).numFrames);
            
            offset_time = calculate_delay_offset(h); % Get offset time
            plot_delayed_signal=0;
            interpolation_factor = 10;
            n=1;
            for n_frame = 1:h.Resource.RcvBuffer(1).numFrames
                for n_tx = 1:length(seq)
                    %% compute time vector for this line
                    t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
                    t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
                    no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
                    
                    % compute the offset in time from center of probe to
                    % center of transmit wave. We do this by finding the
                    % mean between the two center transmit delays for a
                    % even numbered probe, and the center transmit delay
                    % for a odd elemtn probe. We have to calculate the
                    % transmit delays ourselves, since the delays in
                    % Tx.Delay is cropped to only the active elements.
                    trans_delays = calculate_trans_delays(h,channel_data,n_tx);
                    t0_1 = mean(trans_delays(ceil(channel_data.probe.N_elements/2):ceil((channel_data.probe.N_elements+1)/2)));
                    
                    t_out = linspace(t_ini,t_end,no_t);
                    t_in=t_out-offset_time-t0_1;
                    
                    data_in = h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,h.Trans.Connector,n_frame);
                    data(:,:,n_tx,n_frame) = time_shift_data(h,data_in,t_in,t_out,interpolation_factor,channel_data);
                    
                    n=n+1;
                    if plot_delayed_signal
                        %% Point to beamform to (where the scatterer is in the simulation)
                        % Need to change to correct scatter setup in the
                        % Verasonics script, see FI_phase_array_p4.m for
                        % example. This seems to be correct, but the delays
                        % are slighty off for transmit angles > 0 but not
                        % for angles < 0. Strange. Is there somthing wrong
                        % with the Verasonics simulation?? :)
                        
                        
                        %[z_all,x_all] = pol2cart(h.angles,ones(1,length(h.angles))*40e-3);
                        x = channel_data.sequence(n_tx).source.x;
                        y = 0;
                        z = channel_data.sequence(n_tx).source.z;
                        
                        TF=(-1).^(z<channel_data.sequence(n_tx).source.z).*sqrt((channel_data.sequence(n_tx).source.x-x).^2+(channel_data.sequence(n_tx).source.y-y).^2+(channel_data.sequence(n_tx).source.z-z).^2);
                        % add distance from source to origin
                        TF=TF+sign(cos(channel_data.sequence(n_tx).source.azimuth)).*channel_data.sequence(n_tx).source.distance;
                        % receive delay
                        RF=sqrt((channel_data.probe.x-x).^2+(channel_data.probe.y-y).^2+(channel_data.probe.z-z).^2);
                        % total delay
                        delay=(RF+TF)/channel_data.sound_speed;
                        
                        figure(100);plot(h.TX(n_tx).Delay)
                        
                        figure(102); hold off;
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
    end
    
    % Private methods
    methods (Access = private)
        
        % Calculate the offset time from start of transmitted pulse to
        % center of pulse and compensate for the lens correction
        function offset_time = calculate_delay_offset(h)
            % offset calculation
            offset_distance=(h.TW.peak)*h.lambda;   % in [m]
            if strcmp(h.Trans.units,'mm')
                offset_distance=offset_distance+2*h.Trans.lensCorrection*1e-3;
            elseif strcmp(h.Trans.units,'wavelengths')
                offset_distance=offset_distance+2*h.Trans.lensCorrection*h.lambda;
            end
            offset_time=offset_distance/h.c0;   % in [s]
        end
        
        % Generate a USTB probe object from the Verasonics parameters
        function prb = create_probe_object(h)
            if strcmp(h.Trans.name,'L7-4') || strcmp(h.Trans.name,'P4-2v')
                prb=uff.linear_array();
                prb.N=h.Trans.numelements;                  % number of elements
                prb.pitch=h.Trans.spacingMm/1000;           % probe pitch in azimuth [m]
                prb.element_width=h.Trans.elementWidth/1000;   % element width [m]
            else
                error('Sorry, that probe is not supported in USTB yet.');
            end
        end
        
        function trans_delays = calculate_trans_delays(h,channel_data,n_tx)
            %% Stolen from the computeTXDelays ;)
            % Hacked to work for the Verasonics definition of linear_array transmit focus with azimuth = 0
            % Delays are returned in seconds
            angle = 0;
            FocalPt(1) = channel_data.sequence(n_tx).source.x + channel_data.sequence(n_tx).source.z * sin(angle);
            FocalPt(2) = 0.0;
            FocalPt(3) = channel_data.sequence(n_tx).source.z * cos(angle);
            % Compute distance to focal point from each active element.
            X = channel_data.probe.geometry(:,1)' - FocalPt(1);
            D = sqrt(X.*X + FocalPt(3)*FocalPt(3));
            Indices = find(logical(h.TX(n_tx).Apod));
            D = max(D) - D;
            D = D - D(Indices(end));
            
            %figure(101);
            %plot(D); hold on;
            %plot(h.TX(n_tx).Delay*h.lambda)
            
            trans_delays = D/channel_data.sound_speed;
        end
        
        function data_out = time_shift_data(h,data_in,t_in,t_out,interpolation_factor,channel_data)
            % First do a sinc interpolation
            t_in_interp = linspace(t_in(1),t_in(end)+((interpolation_factor-1)/interpolation_factor)*(1/channel_data.sampling_frequency),length(t_in)*interpolation_factor); 
            data_tx_interpolated = interpft(double(data_in),length(t_in)*interpolation_factor);
            %%
            %                     channel = 64;
            %                     figure(99);clf;hold all;
            %                     subplot(211);hold all
            %                     plot(t_in_interp,data_tx_interpolated(:,channel),'Displayname','interpolated');
            %                     plot(t_in,data_tx(:,channel),'Displayname','original');
            %                     subplot(212);hold all
            
            %%
            % read data
            data_out=interp1(t_in_interp,data_tx_interpolated,t_out,'linear',0);
        end
    end
end
