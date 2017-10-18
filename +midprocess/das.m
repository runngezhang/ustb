classdef das < midprocess
    %DAS   Implementation of USTB DAS general beamformer
    %
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %
    %   $Last updated: 2017/09/22$
    
    %% Additional properties
    properties
        dimension = dimension.receive;      % dimension enumeration class that specifies whether the process will run only on transmit, receive, both, or none.
        code = code.mex;                    % code enumeration class that specifies the code to be run (code.matlab, code.mex)
    end
    
    %% constructor
    methods (Access = public)
        function h=das()
            h.name='USTB DAS General Beamformer';
            h.reference= 'www.ustb.no';
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};
            h.version='v1.0.12';
        end
    end
    
    %% go method
    methods
        function beamformed_data=go(h)
            
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
            
            % modulation frequency
            w0=2*pi*h.channel_data.modulation_frequency;
            
            % constants
            sampling_frequency=single(h.channel_data.sampling_frequency);
            initial_time=single(h.channel_data.initial_time);
            modulation_frequency=single(h.channel_data.modulation_frequency);
            
            % calculate transmit apodization according to 10.1109/TUFFC.2015.007183
            h.transmit_apodization.sequence=h.channel_data.sequence;
            h.transmit_apodization.focus=h.scan;
            tx_apodization=single(h.transmit_apodization.data);
            
            % calculate receive apodization
            h.receive_apodization.probe=h.channel_data.probe;
            h.receive_apodization.focus=h.scan;
            rx_apodization=single(h.receive_apodization.data);
            rx_propagation_distance=(h.receive_apodization.propagation_distance);
            
            % calculate receive delay
            xm=bsxfun(@minus,h.channel_data.probe.x.',h.scan.x);
            ym=bsxfun(@minus,h.channel_data.probe.y.',h.scan.y);
            zm=bsxfun(@minus,h.channel_data.probe.z.',h.scan.z);
            receive_delay=single(sqrt(xm.^2+ym.^2+zm.^2)/h.channel_data.sound_speed);
            
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
            transmit_delay = single(transmit_delay./h.channel_data.sound_speed);
            
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
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %           MEX
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case code.mex
                        aux_data=mex.das_c(data,...
                                           sampling_frequency,...
                                           initial_time,...
                                           tx_apodization,...
                                           rx_apodization,...
                                           transmit_delay,...
                                           receive_delay,...
                                           modulation_frequency,...
                                           int32(h.dimension));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %        MATLAB
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case code.matlab %#ok<PROP>
                        % workbar
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
                                        if mod(n,round(N/100))==0
                                            tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB MATLAB');
                                        end
                                        
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
                                            case dimension.none %#ok<*PROP>
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
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %        MATLAB GPU
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case code.matlab_gpu
                        
                        % set auxiliary data to be a gpuArray
                        
                        if (h.code==code.matlab_gpu)
                            gpu_info=gpuDevice();
                            aux_data_memory=prod([size(aux_data) 8]) + prod([N_pixels N_waves N_frames 8]);
                            if (aux_data_memory*1.3>gpu_info.AvailableMemory)
                                warning(sprintf('DAS GPU: The beamformed data %0.2f MB does not fit in the available GPU memory %0.2f MB. GPU computation will performed sequentially.',aux_data_memory/2^20, gpu_info.AvailableMemory/2^20));
                            else
                                fprintf(1, 'DAS GPU: The beamformed data %0.2f MB fits in the available GPU memory %0.2f MB.',aux_data_memory/2^20, gpu_info.AvailableMemory/2^20);
                                aux_data=gpuArray(aux_data);
                            end
                        end
                        
                        % workbar
                        tools.workbar();
                        N=N_waves*N_channels;
                        
                        % convert structures to gpuArray
                        time_gpu= gpuArray(h.channel_data.time);
                        w0_gpu= gpuArray(w0);

                        % receive loop
                        for n_rx=1:N_channels
                            if any(rx_apodization(:,n_rx))

                                % move chunk of channel data to array
                                data_gpu = gpuArray(data(:,n_rx,:,:));

                                % transmit loop
                                for n_wave=1:N_waves
                                    if any(tx_apodization(:,n_wave))
                                        
                                        % progress bar
                                        n=(n_rx-1)*N_waves+n_wave;
                                        if mod(n,round(N/100))==0
                                            tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB GPU');
                                        end
                                        
                                        apodization_gpu= gpuArray(rx_apodization(:,n_rx).*tx_apodization(:,n_wave));
                                        delay_gpu= gpuArray(receive_delay(:,n_rx) + transmit_delay(:,n_wave));
                                        
                                        % apply phase correction factor to IQ data
                                        if(w0>eps) 
                                            apodization_gpu = exp(1i.*w0_gpu*delay_gpu).*apodization_gpu;
                                        end
                                        
                                        % beamformed signal
                                        temp = bsxfun(@times,apodization_gpu,interp1(time_gpu,data_gpu(:,1,n_wave,:),delay_gpu,'linear',0));
                                        
                                        % set into auxiliary data
                                        if ~isa(aux_data,'gpuArray') temp=gather(temp); end
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
                        
                        % assign phase according to 2 times the receive propagation distance
                        if(w0>eps)
                            aux_data=bsxfun(@times,aux_data,exp(-1i*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
                        end
                        
                        % gather & clear GPU memory
                        if(h.code == code.matlab_gpu || h.code == code.matlab_gpu_frameloop)
                            aux_data=gather(aux_data);
                            gpuDevice(1);
                        end
                        
                        tools.workbar(1);
                        
                    % =====================================
                    %
                    %         MATLAB GPU FRAMELOOP
                    %
                    % =====================================
                    case code.matlab_gpu_frameloop
                        gpu_info=gpuDevice();
                        fprintf(1, '=========== DAS GPU FRAMELOOP ===========\n');
                        fprintf(1, 'Selected GPU: %s\n', gpu_info.Name);
                        fprintf(1, 'Available memory: %0.2f MB\n', gpu_info.AvailableMemory*1e-6);
                        
                        % transfer data to the GPU
                        time_gpu        = gpuArray(single(h.channel_data.time));
                        w0_gpu          = gpuArray(single(w0));
                        c0_gpu          = gpuArray(single(h.channel_data.sound_speed));
                        rx_apod_gpu     = gpuArray(rx_apodization);
                        tx_apod_gpu     = gpuArray(tx_apodization);
                        rx_delay_gpu    = gpuArray(receive_delay);
                        tx_delay_gpu    = gpuArray(transmit_delay);
                        rx_distance_gpu = gpuArray(single(rx_propagation_distance));
                        phase_term_gpu  = exp(-1j*w0_gpu*2*rx_distance_gpu/c0_gpu);
                        
                        if N_waves == 1
                            apod_gpu  = bsxfun(@times,rx_apod_gpu,tx_apod_gpu);
                            delay_gpu = bsxfun(@plus, rx_delay_gpu, tx_delay_gpu);
                        end
                        
                        % workbar
                        tools.workbar();
                        
                        % frame loop
                        for n_frame=1:N_frames
                            % update workbar
                            tools.workbar(n_frame/N_frames, sprintf('%s (%s)',h.name, h.version),'USTB');
                            
                            % transfer channel data to device
                            ch_data = gpuArray(single(data(:,:,:,n_frame)));
                            
                            % beamformed data buffer
                            switch h.dimension
                                case dimension.transmit
                                    bf_data = complex(zeros([N_pixels, N_channels], 'single', 'gpuArray'));
                                case dimension.both
                                    bf_data = complex(zeros([N_pixels, 1, N_waves], 'single', 'gpuArray'));
                            end
                            
                            % wave loop
                            for n_wave=1:N_waves
                                
                                if N_waves > 1
                                    apod_gpu  = bsxfun(@times,rx_apod_gpu,tx_apod_gpu(:,n_wave));
                                    delay_gpu = bsxfun(@plus, rx_delay_gpu, tx_delay_gpu(:,n_wave));
                                end
                                
                                % apply phase correction factor to IQ data
                                if(w0>eps)
                                    apod_gpu = exp(1i.*w0_gpu*delay_gpu).*apod_gpu;
                                end
                                
                                pre_bf_data = complex(zeros([N_pixels, N_channels], 'single', 'gpuArray'));
                                
                                % channel loop
                                for n_rx=1:N_channels
                                    pre_bf_data(:,n_rx) = interp1(time_gpu, ch_data(:,n_rx, n_wave, :), delay_gpu(:,n_rx), 'linear',0);
                                end % end channel loop
                                
                                % apply apodization
                                pre_bf_data = apod_gpu .* pre_bf_data;
                                
                                % assign phase according to 2 times the receive propagation distance
                                if(w0_gpu>eps)
                                    pre_bf_data = pre_bf_data .* phase_term_gpu;
                                end
                                
                                switch h.dimension
                                    case dimension.none
                                        aux_data(:,n_rx,n_wave,n_frame) = gather(pre_bf_data);
                                    case dimension.receive
                                        aux_data(:,1,n_wave,n_frame) = gather(sum(pre_bf_data, 2));
                                    case dimension.transmit
                                        bf_data = bf_data + pre_bf_data;
                                    case dimension.both
                                        bf_data(:,1,N_waves) = sum(pre_bf_data, 2);
                                end           
                            end % end wave loop
                            
                            switch h.dimension
                                case dimension.transmit
                                    aux_data(:,:,1,n_frame) = gather(sum(bf_data, 3));
                                case dimension.both
                                    aux_data(:,1,1,n_frame) = gather(sum(bf_data, 3));
                            end
                        end % end frame loop
                        
                        tools.workbar(1);
                        
                    % =======================================
                    %
                    %       MATLAB GPU FRAMELOOP CHUNK
                    %
                    % =======================================
                    case code.matlab_gpu_frameloop_chunk
                        
                        % check available GPU memory and set the chunk size accordingly
                        gpu_info=gpuDevice();
                        rx_delay_size = numel(receive_delay)  * 4;
                        tx_delay_size = numel(receive_delay)  * 4;
                        rx_apod_size  = numel(rx_apodization) * 4;
                        tx_apod_size  = numel(tx_apodization) * 4;
                        ch_data_size  = numel(size(data(:,:,:,1))) * 8;
                        
                        availableMemory = @(x)  gpu_info.AvailableMemory - ...
                                                rx_delay_size - ...
                                                tx_delay_size - ...
                                                rx_apod_size - ...
                                                tx_apod_size - ...
                                                ch_data_size*x;
                        chunk_size = min(   [fzero(@(x) floor(availableMemory(x)/(N_pixels*8)) - 2*x, 1); ...
                                            floor(intmax('int32')/(N_pixels)); ...
                                            N_frames]);
                        
                        
                        fprintf(1, '=========== DAS GPU FRAMELOOP CHUNK ===========\n');
                        fprintf(1, 'Selected GPU: %s\n', gpu_info.Name);
                        fprintf(1, 'Available memory: %0.2f MB\n', gpu_info.AvailableMemory*1e-6);
                        fprintf(1, 'GPU can accomodate %d frames at a time\n', chunk_size);
                        
                        % workbar
                        tools.workbar();
                        
                        % transfer data to the GPU
                        time_gpu        = gpuArray(single(h.channel_data.time));
                        w0_gpu          = gpuArray(single(w0));
                        c0_gpu          = gpuArray(single(h.channel_data.sound_speed));
                        rx_apod_gpu     = gpuArray(rx_apodization);
                        tx_apod_gpu     = gpuArray(tx_apodization);
                        rx_delay_gpu    = gpuArray(receive_delay);
                        tx_delay_gpu    = gpuArray(transmit_delay);
                        rx_distance_gpu = gpuArray(single(rx_propagation_distance));
                        phase_term_gpu  = exp(-1j*w0_gpu*2*rx_distance_gpu/c0_gpu);
                        
                        % frame loop
                        for n_frame=1:chunk_size:N_frames
                            % channel data chunk
                            chunk_start = n_frame;
                            
                            if mod(N_frames-n_frame, chunk_size) == 0
                                chunk_end = n_frame+chunk_size-1;
                            else
                                chunk_end = N_frames;
                            end
                            
                            ch_data = gpuArray(single(data(:,:,:,chunk_start:chunk_end)));
                            
                            % beamformed data buffer
                            bf_data = complex(zeros([N_pixels, 1, 1, chunk_size], 'single', 'gpuArray'));
                            
                            % transmit loop
                            for n_wave=1:N_waves
                                tools.workbar((n_wave + (n_frame-1))/N_waves/N_frames, sprintf('%s (%s)',h.name, h.version),'USTB GPU chunk frameloop');
                                if any(tx_apod_gpu(:,n_wave))
                                    
                                    apod_gpu  = bsxfun(@times,rx_apod_gpu,tx_apod_gpu(:,n_wave));
                                    delay_gpu = bsxfun(@plus, rx_delay_gpu, tx_delay_gpu(:,n_wave));
                                    
                                    % apply phase correction factor to IQ data
                                    if(w0>eps)
                                        apod_gpu = exp(1i.*w0_gpu*delay_gpu).*apod_gpu;
                                    end
                                    
                                    % channel loop
                                    for n_rx=1:N_channels
                                        if any(rx_apod_gpu(:,n_rx))
                                            pre_bf_data = interp1(time_gpu, ch_data(:,n_rx, n_wave, :), delay_gpu(:,n_rx), 'linear',0);
                                            
                                            % apply apodization
                                            pre_bf_data = bsxfun(@times, apod_gpu(:,n_rx), pre_bf_data);
                                            
                                            % assign phase according to 2 times the receive propagation distance
                                            if(w0_gpu>eps)
                                                pre_bf_data = bsxfun(@times, pre_bf_data, phase_term_gpu);
                                            end
                                            
                                            switch h.dimension
                                                case dimension.none
                                                    aux_data(:,n_rx,n_wave,chunk_start:chunk_end) = gather(pre_bf_data);
                                                case dimension.receive
                                                    bf_data = bf_data + pre_bf_data;
                                                case dimension.transmit
                                                    aux_data(:,:,n_wave,chunk_start:chunk_end)
                                                case dimension.both
                                                    bf_data = bf_data + pre_bf_data;
                                            end
                                        end
                                    end % end channel loop
                                    
                                    if h.dimension == dimension.receive
                                        aux_data(:,1,n_wave,chunk_start:chunk_end) = gather(bf_data);
                                    end
                                end
                            end % end waves loop
                            if h.dimension == dimension.both
                                aux_data(:,1,1,chunk_start:chunk_end) = gather(bf_data);
                            end
                        end % end frame loop
                        
                        tools.workbar(1);
                  
                    otherwise
                        error('Unknown code implementation requested');
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
