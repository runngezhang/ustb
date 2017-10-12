classdef gpu_das_matlab < midprocess
    %   GPU_DAS_MATLAB  Matlab implementation of the Delay-and-Sum general
    %                   beamformed with gpu
    %
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %            Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %            Stefano Fiorentini     <stefano.fiorentini@ntnu.no>
    %
    %   $Last updated: 2017/10/5$
    
    %% constructor
    methods (Access = public)
        function h=gpu_das_matlab()
            h.name='USTB GPU enabled, general DAS Beamformer MATLAB';
            h.reference= 'www.ustb.no';
            h.implemented_by={'Stefano Fiorentini <stefano.fiorentini@ntnu.no>', 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};
            h.version='v1.0.0';
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
            
            % modulation frequency
            w0 = 2*pi*gpuArray(h.channel_data.modulation_frequency);
            c0 = gpuArray(h.channel_data.sound_speed);
            
            [channel_data_time, channel_idx_in] = ndgrid(gpuArray(h.channel_data.time), gpuArray(1:h.channel_data.N_channels));
            channel_idx_out = repmat(gpuArray(1:h.channel_data.N_channels), [h.scan.N_pixels, 1]);
            
            % precalculate receive apodization
            h.receive_apodization.probe=h.channel_data.probe;
            h.receive_apodization.focus=h.scan(1);
            rx_apo=gpuArray(h.receive_apodization.data);
            rx_propagation_distance=gpuArray(h.receive_apodization.propagation_distance);
            
            % precalculate transmit apodization according to 10.1109/TUFFC.2015.007183
            h.transmit_apodization.sequence=h.channel_data.sequence;
            h.transmit_apodization.focus=h.scan(1);
            tx_apodization=gpuArray(h.transmit_apodization.data);
            
            % precalculate receive delay
            probe_x = gpuArray(h.channel_data.probe.x).';
            probe_y = gpuArray(h.channel_data.probe.y).';
            probe_z = gpuArray(h.channel_data.probe.z).';
            
            xm=bsxfun(@minus,probe_x,gpuArray(h.scan(1).x));
            ym=bsxfun(@minus,probe_y,gpuArray(h.scan(1).y));
            zm=bsxfun(@minus,probe_z,gpuArray(h.scan(1).z));
            RF=sqrt(xm.^2+ym.^2+zm.^2);
            
            % precalculating hilbert (if needed)
            data=h.channel_data.data;
            if ~(w0>eps)
                data=reshape(hilbert(h.channel_data.data(:,:)),size(h.channel_data.data));
            end
            
            % create beamformed data class
            h.beamformed_data=uff.beamformed_data();
            h.beamformed_data.scan=h.scan;
            N_pixels = 0; for n=1:length(h.scan) N_pixels = max([N_pixels h.scan(n).N_pixels]); end
            % out_data.sequence=h.channel_data.sequence; % not included by default
            
            % auxiliary data
            bf_data=zeros(N_pixels,1,numel(h.channel_data.sequence),h.channel_data.N_frames);
            
            % wave loop
            tools.workbar();
            N=numel(h.channel_data.sequence)*h.channel_data.N_frames;
            for n_wave=1:numel(h.channel_data.sequence)
                
                % support multiple or single scans with the same code
                if numel(h.scan)==1
                    current_scan=h.scan;
                else
                    current_scan=h.scan(n_wave);
                end         
                current_scan_x = gpuArray(current_scan.x);
                current_scan_y = gpuArray(current_scan.y);
                current_scan_z = gpuArray(current_scan.z);
                
                % calculate receive apodization for multiple scan
                if numel(h.scan)>1
                    h.receive_apodization.focus=current_scan;
                    rx_apo=gpuArray(h.receive_apodization.data);
                    rx_propagation_distance=gpuArray(h.receive_apodization.propagation_distance);
                end
                
                % calculate transmit apodization for multiple scan
                if numel(h.scan)>1
                    h.transmit_apodization.sequence=h.channel_data.sequence(n_wave);
                    h.transmit_apodization.focus=current_scan;
                    tx_apo=gpuArray(h.transmit_apodization.data);
                else
                    tx_apo=tx_apodization(:,n_wave);
                end   
                
                % transmit delay
                source_x        = gpuArray(h.channel_data.sequence.source.x);
                source_y        = gpuArray(h.channel_data.sequence.source.y);
                source_z        = gpuArray(h.channel_data.sequence.source.z);
                source_azimuth  = gpuArray(h.channel_data.sequence.source.azimuth);
                source_elevation= gpuArray(h.channel_data.sequence.source.elevation);
                source_distance = gpuArray(h.channel_data.sequence.source.distance);
                
                if ~isinf(source_distance)
                    % point sources
                    TF=(-1).^(current_scan_z<source_z).*sqrt((source_x-current_scan_x).^2+(source_y-current_scan_y).^2+(source_z-current_scan_z).^2);
                    % add distance from source to origin
                    %OLD VERSION with rounding problem: TF=TF+sign(cos(h.channel_data.sequence(n_wave).source.azimuth)).*h.channel_data.sequence(n_wave).source.distance;
                    if (source_z<-1e-3)
                        TF=TF-source_distance;
                    else
                        TF=TF+source_distance;
                    end
                else
                    % plane waves
                    TF=current_scan_z*cos(source_azimuth)*cos(source_elevation)+current_scan_x*sin(source_azimuth)*cos(source_elevation)+current_scan_y*sin(source_elevation);
                end
                
                % calculate receive delay for multiple scan
                if numel(h.scan)>1
                    xm=bsxfun(@minus,probe_x,current_scan_x);
                    ym=bsxfun(@minus,probe_y,current_scan_y);
                    zm=bsxfun(@minus,probe_z,current_scan_z);
                    RF=sqrt(xm.^2+ym.^2+zm.^2);
                end
                

                
                % total delay
                delay=bsxfun(@plus, RF, TF)/c0;
                
                for n_frame=1:h.channel_data.N_frames
                    % progress bar
                    tools.workbar((n_frame + (n_wave-1)*h.channel_data.N_frames)/N, sprintf('%s (%s)',h.name, h.version),'USTB');
                    
                    % phase correction factor
                    if(w0>eps)
                        phase_shift=exp(1i.*w0*delay);
                    else
                        phase_shift=1;
                    end
                    
                    % beamformed signal
                    ch_data = gpuArray(data(:,:,n_wave,n_frame));
                    pre_bf_data = interpn(  channel_data_time, channel_idx_in, ...
                                            ch_data, ...
                                            delay, channel_idx_out, 'linear');
                                        
                    bf_data(:,1,n_wave,n_frame) = gather(sum(tx_apo.*rx_apo.*phase_shift.*pre_bf_data, 2, 'omitnan').*exp(-1j*w0*2*rx_propagation_distance/c0));
                end
            end
            h.beamformed_data.data=bf_data;
            tools.workbar(1);
            
            % pass reference
            beamformed_data = h.beamformed_data;
            
            % update hash
            h.save_hash();
        end
    end
end
