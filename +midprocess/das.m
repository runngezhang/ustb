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

            % clear GPU memory
            if(h.code == code.matlab_gpu || h.code == code.matlab_gpu_frameloop) gpuDevice(1); end
            
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
            if(h.code==code.matlab_gpu || h.code==code.matlab_gpu_frameloop ) 
                rx_propagation_distance=gpuArray(rx_propagation_distance); 
            end
            
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
            % set auxiliary data to be a gpuArray
            if (h.code==code.matlab_gpu)
                gpu_info=gpuDevice();
                aux_data_memory=prod([size(aux_data) 8]) + prod([N_pixels N_waves N_frames 8]);
                if (aux_data_memory*1.3>gpu_info.AvailableMemory)
                    warning(sprintf('DAS GPU: The beamformed data %0.2f MB does not fit in the available GPU memory %0.2f MB. GPU computation will performed sequentially.',aux_data_memory/2^20, gpu_info.AvailableMemory/2^20));
                else
                    disp(sprintf('DAS GPU: The beamformed data %0.2f MB fits in the available GPU memory %0.2f MB.',aux_data_memory/2^20, gpu_info.AvailableMemory/2^20));
                    aux_data=gpuArray(aux_data);
                end
            end
            % set auxiliary data to be a gpuArray
            if (h.code==code.matlab_gpu_frameloop)
                gpu_info=gpuDevice();
                aux_data_memory=prod([size(aux_data) 8]) + prod([N_pixels N_channels 8]);
                if (aux_data_memory*1.3>gpu_info.AvailableMemory)
                    warning(sprintf('DAS GPU: The beamformed data %0.2f MB does not fit in the available GPU memory %0.2f MB. GPU computation will performed sequentially.',aux_data_memory/2^20, gpu_info.AvailableMemory/2^20));
                else
                    disp(sprintf('DAS GPU: The beamformed data %0.2f MB fits in the available GPU memory %0.2f MB.',aux_data_memory/2^20, gpu_info.AvailableMemory/2^20));
                    aux_data=gpuArray(aux_data);
                end
            end
            
            % delay & sum
            if any(data(:)>0) % only process if any data > 0
                
                switch h.code
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% MEX
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
                    %%% MATLAB
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case code.matlab
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
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% MATLAB GPU
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case code.matlab_gpu
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
                                        
                                        % clear gpu variables
                                        clear temp delay_gpu apodization_gpu;
                                    end
                                end
                            end
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% MATLAB GPU FRAMELOOP
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case code.matlab_gpu_frameloop
                        % workbar
                        tools.workbar();
                        N=N_waves*N_channels;

                        % convert structures to gpuArray
                        time_gpu= gpuArray(h.channel_data.time);
                        w0_gpu= gpuArray(w0);

                        % transmit loop
                        for n_wave=1:N_waves
                            if any(tx_apodization(:,n_wave))
                                
                                apodization_gpu= gpuArray(bsxfun(@times,rx_apodization,tx_apodization(:,n_wave)));
                                delay_gpu= gpuArray(bsxfun(@plus, receive_delay, transmit_delay(:,n_wave)));
                                
                                % apply phase correction factor to IQ data
                                if(w0>eps) 
                                    apodization_gpu = exp(1i.*w0_gpu*delay_gpu).*apodization_gpu;
                                end
                                
                                % frame loop
                                for n_frame=1:N_frames
                                    % progress bar
                                    tools.workbar((n_frame + (n_wave-1)*N_frames)/N_frames/N_waves, sprintf('%s (%s)',h.name, h.version),'USTB GPU frameloop');
                                
                                    % channel data
                                    data_gpu= gpuArray(single(data(:,:,n_wave,n_frame)));

                                    % temporal data
                                    temp = zeros([h.scan.N_pixels, h.channel_data.N_channels], 'single', 'gpuArray');
    
                                    % channel loop
                                    for n_rx=1:N_channels
                                        if any(rx_apodization(:,n_rx))
                                            temp(:, n_rx) = interp1(time_gpu,data_gpu(:,n_rx),delay_gpu(:,n_rx),'linear',0);
                                        end
                                    end
                                    
                                    % apply apodization
                                    temp = bsxfun(@times,apodization_gpu,temp);
                                            
                                    % set into auxiliary data
                                    if ~isa(aux_data,'gpuArray') temp=gather(temp); end                                    
                                    switch h.dimension
                                        case dimension.none
                                            aux_data(:,:,n_wave,n_frame)=temp;
                                        case dimension.receive
                                            aux_data(:,1,n_wave,n_frame)=aux_data(:,1,n_wave,n_frame)+sum(temp,2);
                                        case dimension.transmit
                                            aux_data(:,:,1,n_frame)=aux_data(:,n_rx,1,n_frame)+temp;
                                        case dimension.both
                                            aux_data(:,1,1,n_frame)=aux_data(:,1,1,n_frame)+sum(temp,2);
                                    end
                                    
                                    % clear gpu variables
                                    clear temp data_gpu;
                                end
                            end
                        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    otherwise
                        error('Unknown code implementation requested');
                end
                tools.workbar(1);

                % assign phase according to 2 times the receive propagation distance
                if(w0>eps)
                    aux_data=bsxfun(@times,aux_data,exp(-1i*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
                end
                
                % gather & clear GPU memory
                if(h.code == code.matlab_gpu || h.code == code.matlab_gpu_frameloop) 
                    aux_data=gather(aux_data); 
                    gpuDevice(1); 
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
