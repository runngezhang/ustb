classdef das_matlab_gpu_nvidia < midprocess
    %DAS   Implementation of USTB DAS general beamformer
    %
    %   authors: Stefano Fiorentini     <stefano.fiorentini@ntnu.no>
    %            Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %            Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %
    %   $Last updated: 2017/10/5$
    
    %% Additional properties
    properties
        dimension = dimension.receive;      % dimension enumeration class that specifies whether the process will run only on transmit, receive, both, or none.
    end
    
    %% constructor
    methods (Access = public)
        function h=das()
            h.name='USTB DAS MATLAB GPU Nvidia enabled';
            h.reference= 'www.ustb.no';
            h.implemented_by={'Stefano Fiorentini <stefano.fiorentini@ntnu.no>', 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};
            h.version='v1.0.1';
        end
    end
    
    %% go method
    methods
        function beamformed_data=go(h)
            % clear GPU memory
            gpuDevice(1);
            
            % check if we can skip calculation
            if h.check_hash()
                beamformed_data= h.beamformed_data;
                return;
            end
            
            % short names
            N_pixels = h.scan.N_pixels;
            N_channels = h.channel_data.N_channels;
            N_waves = h.channel_data.N_waves;
            N_frames = h.channel_data.N_frames;
            
            % defining some variables gpuArrays
            %w0=2*pi*h.channel_data.modulation_frequency;
            w0 = 2*pi*gpuArray(h.channel_data.modulation_frequency);
            c0 = gpuArray(h.channel_data.sound_speed);
            [channel_data_time, channel_idx_in] = ndgrid(gpuArray(h.channel_data.time), gpuArray(1:h.channel_data.N_channels));
            channel_idx_out = repmat(gpuArray(1:h.channel_data.N_channels), [h.scan.N_pixels, 1]);
            
            % constants
            sampling_frequency=single(h.channel_data.sampling_frequency);
            initial_time=single(h.channel_data.initial_time);
            modulation_frequency=single(h.channel_data.modulation_frequency);
            
            % calculate transmit apodization according to 10.1109/TUFFC.2015.007183
            h.transmit_apodization.sequence=h.channel_data.sequence;
            h.transmit_apodization.focus=h.scan;
            tx_apodization=gpuArray(h.transmit_apodization.data);
            
            % calculate receive apodization
            h.receive_apodization.probe=h.channel_data.probe;
            h.receive_apodization.focus=h.scan;
            rx_apodization=gpuArray(h.receive_apodization.data);
            rx_propagation_distance=gpuArray(h.receive_apodization.propagation_distance);
            
            % calculate receive delay
            xm=bsxfun(@minus,h.channel_data.probe.x.',h.scan.x);
            ym=bsxfun(@minus,h.channel_data.probe.y.',h.scan.y);
            zm=bsxfun(@minus,h.channel_data.probe.z.',h.scan.z);
            receive_delay=gpuArray(sqrt(xm.^2+ym.^2+zm.^2)/h.channel_data.sound_speed);
            
            % calculate transmit delay
            transmit_delay=zeros(N_pixels,N_waves);
            for n_wave=1:numel(h.channel_data.sequence)
                if ~isinf(h.channel_data.sequence(n_wave).source.distance)
                    % point sources
                    transmit_delay(:,n_wave)=(-1).^(h.scan.z<h.channel_data.sequence(n_wave).source.z).*sqrt((h.channel_data.sequence(n_wave).source.x-h.scan.x).^2+(h.channel_data.sequence(n_wave).source.y-h.scan.y).^2+(h.channel_data.sequence(n_wave).source.z-h.scan.z).^2);
                    
                    % add distance from source to origin
                    if (h.channel_data.sequence(n_wave).source.z<-1e-3)
                        transmit_delay(:,n_wave)=transmit_delay(:,n_wave)-h.channel_data.sequence(n_wave).source.distance;
                    else
                        transmit_delay(:,n_wave)=transmit_delay(:,n_wave)+h.channel_data.sequence(n_wave).source.distance;
                    end
                else
                    % plane waves
                    transmit_delay(:,n_wave)=h.scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
                end
            end
            transmit_delay = gpuArray(transmit_delay./h.channel_data.sound_speed);
            
            % precalculating hilbert (if needed)
            tools.check_memory(prod([size(h.channel_data.data) 8]));
            data=single(h.channel_data.data);
            if ~(w0>eps)
                data=single(reshape(hilbert(h.channel_data.data(:,:)),size(h.channel_data.data)));
            end
            
            % create beamformed data class
            h.beamformed_data=uff.beamformed_data();
            h.beamformed_data.scan=h.scan;
            
            % auxiliary data
            switch h.dimension
                case dimension.none
                    tools.check_memory(prod([N_pixels N_channels N_waves N_frames 8]));
                    aux_data=zeros(N_pixels,N_channels,N_waves,N_frames);
                case dimension.receive
                    tools.check_memory(prod([N_pixels N_waves N_frames 8]));
                    aux_data=zeros(N_pixels,1,N_waves,N_frames);
                case dimension.transmit
                    tools.check_memory(prod([N_pixels N_channels N_frames 8]));
                    aux_data=zeros(N_pixels,N_channels,1,N_frames);
                case dimension.both
                    tools.check_memory(prod([N_pixels N_frames 8]));
                    aux_data=zeros(N_pixels,1,1,N_frames);
            end
            
            % delay & sum
            if any(data(:)>0) % only process if any data > 0
                
                switch h.code
%                     case code.mex
%                         aux_data=mex.das_c(data,...
%                                            sampling_frequency,...
%                                            initial_time,...
%                                            tx_apodization,...
%                                            rx_apodization,...
%                                            transmit_delay,...
%                                            receive_delay,...
%                                            modulation_frequency,...
%                                            int32(h.dimension));
                    case code.matlab
                        % apply phase shift to apodization
                        tools.workbar();
                        N=N_waves*N_channels;
                        
                        % transmit loop
                        for n_wave=1:N_waves
                            if any(tx_apodization(:,n_wave))
                            
                                % receive loop
                                for n_rx=1:N_channels
                                    if any(rx_apodization(:,n_rx))
                                        
                                        % progress bar
                                        n=(n_wave-1)*N_channels+n_rx;
                                        %if mod(n,round(N/100))==0
                                            tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                                        %end
                                        
                                        apodization= rx_apodization(:,n_rx).*tx_apodization(:,n_wave);
                                        delay= receive_delay(:,n_rx) + transmit_delay(:,n_wave);
                                        
                                        % beamformed signal
                                        temp = bsxfun(@times,apodization,interp1(h.channel_data.time,data(:,n_rx,n_wave,:),delay,'linear',0));
                                        
                                        % apply phase correction factor to IQ data
                                        if(w0>eps) 
                                            temp = bsxfun(@times,exp(1i.*w0*delay),temp);
                                        end
                                        
                                        % set into auxiliary data
                                        switch h.dimension
                                            case dimension.none
                                                aux_data(:,n_rx,n_wave,:)=temp;
                                            case dimension.receive
                                                aux_data(:,1,n_wave,:)=aux_data(:,1,n_wave,:)+temp;
                                            case dimension.transmit
                                                aux_data(:,n_rx,1,:)=aux_data(:,n_rx,1,:)+temp;
                                            case dimension.both
                                                aux_data(:,1,1,:)=aux_data(:,1,1,:)+temp;
                                        end
                                    end
                                end
                            end
                        end
                    otherwise
                        error('Unknown code implementation requested');
                end
                tools.workbar(1);
                
                % assign phase according to 2 times the receive propagation distance
                if(w0>eps)
                    aux_data=bsxfun(@times,aux_data,exp(-1i*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
                end
            end
            
            % copy data to object
            h.beamformed_data.data = aux_data;
            
            % pass a reference
            beamformed_data = h.beamformed_data;
            
            % update hash
            h.save_hash();
        end
    end
end
