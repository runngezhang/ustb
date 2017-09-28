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
            h.version='v1.0.11';
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
            
            % multiple scans
            if numel(h.scan)>1
                assert(numel(h.scan)==N_waves,'If multiples scans are provided the number of scans must match the number of waves in the sequence.');
                
                % calculate transmit apodization
                tx_apodization=zeros(h.scan(1).N_pixels*N_waves,N_waves);
                stacked_scan = uff.scan();
                for ntx=1:numel(h.scan)
                    assert(h.scan(1).N_pixels==h.scan(ntx).N_pixels,'If multiple scans are provided then all must have the same number of pixels.');
                    stacked_scan.xyz = [stacked_scan.xyz; h.scan(ntx).xyz];
                    tx_apodization((1:N_pixels)+(ntx-1)*N_pixels,ntx)=ones(N_pixels,1);
                end
                h.scan = stacked_scan;
                
                % update number of pixels
                N_pixels =  h.scan.N_pixels;
                
                % change dimension
                if h.dimension == dimension.none
                    h.dimension = dimension.transmit;
                elseif h.dimension == dimension.receive
                    h.dimension = dimension.both;
                end
            else
                % precalculate transmit apodization according to 10.1109/TUFFC.2015.007183
                h.transmit_apodization.sequence=h.channel_data.sequence;
                h.transmit_apodization.focus=h.scan;
                tx_apodization=h.transmit_apodization.data;
            end
            
            % precalculate receive apodization
            h.receive_apodization.probe=h.channel_data.probe;
            h.receive_apodization.focus=h.scan;
            rx_apodization=h.receive_apodization.data;
            rx_propagation_distance=h.receive_apodization.propagation_distance;
            
            % precalculate apodization
            tools.check_memory(prod([size(rx_apodization) size(tx_apodization,2) 8]));
            apodization = single(bsxfun(@times,rx_apodization,permute(tx_apodization,[1 3 2])));
            
            % precalculate receive delay
            xm=bsxfun(@minus,h.channel_data.probe.x.',h.scan.x);
            ym=bsxfun(@minus,h.channel_data.probe.y.',h.scan.y);
            zm=bsxfun(@minus,h.channel_data.probe.z.',h.scan.z);
            RF=sqrt(xm.^2+ym.^2+zm.^2);
            
            % precalculate transmit delay
            TF=zeros(N_pixels,1,N_waves);
            for n_wave=1:numel(h.channel_data.sequence)
                if ~isinf(h.channel_data.sequence(n_wave).source.distance)
                    % point sources
                    TF(:,1,n_wave)=(-1).^(h.scan.z<h.channel_data.sequence(n_wave).source.z).*sqrt((h.channel_data.sequence(n_wave).source.x-h.scan.x).^2+(h.channel_data.sequence(n_wave).source.y-h.scan.y).^2+(h.channel_data.sequence(n_wave).source.z-h.scan.z).^2);
                    % add distance from source to origin
                    if (h.channel_data.sequence(n_wave).source.z<-1e-3)
                        TF(:,1,n_wave)=TF(:,1,n_wave)-h.channel_data.sequence(n_wave).source.distance;
                    else
                        TF(:,1,n_wave)=TF(:,1,n_wave)+h.channel_data.sequence(n_wave).source.distance;
                    end
                else
                    % plane waves
                    TF(:,1,n_wave)=h.scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
                end
            end
            
            % precalculate delay
            tools.check_memory(prod([size(RF) size(TF,2) 8]));
            delay=single(bsxfun(@plus,RF,TF)./h.channel_data.sound_speed);
            
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
                    case code.mex
                        aux_data=mex.das_c(data,sampling_frequency,initial_time,delay,apodization,modulation_frequency,int32(h.dimension));
                    case code.matlab
                        % apply phase shift to apodization
                        if(w0>eps)
                            apodization=exp(1i.*w0*delay).*apodization;
                        end
                        tools.workbar();
                        N=N_waves*N_channels;
                        
                        for n_wave=1:N_waves;
                            for n_rx=1:N_channels
                                % progress bar
                                n=(n_wave-1)*N_channels+n_rx;
                                if mod(n,round(N/100))==1
                                    tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                                end
                                if any(apodization(:,n_rx,n_wave)>0)
                                    % beamformed signal
                                    temp = bsxfun(@times,apodization(:,n_rx,n_wave),interp1(h.channel_data.time,data(:,n_rx,n_wave,:),delay(:,n_rx,n_wave),'linear',0));

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
