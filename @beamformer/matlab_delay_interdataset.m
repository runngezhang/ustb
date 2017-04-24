        function inter_dataset=matlab_delay_interdataset(h,postprocess)
            % MATLAB implementation of the general beamformer
            %
            %   Update: Added modification to be able to do adaptive_beamforming
            
            %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
            %            Ole Marius Hoel Rindal (olemarius@olemarius.net)
            %
            %   $Date: 2017/03/14$
            %          2017/04/07   % Added adaptive beamforming.
            
            
            name='USTB General Beamformer MATLAB'; % beamformer name
            version='v1.0.5';                      % & version

            % modulation frequency
            w0=2*pi*h.channel_data.modulation_frequency;
            
            %% beamforming
            counter = 1;
            total_count = (numel(h.channel_data.sequence)*h.channel_data.N_frames*h.channel_data.N_elements);
            tools.workbar(1/total_count,sprintf('%s (%s)',name,version),'USTB');
            
            % wave loop
            for n_wave=1:numel(h.channel_data.sequence)
                
                % support multiple or single scans with the same code
                if numel(h.scan)==1
                    current_scan=h.scan;
                else
                    current_scan=h.scan(n_wave);
                end
                
                % precalculate receive apodization
                h.receive_apodization.probe=h.channel_data.probe;
                h.receive_apodization.scan=current_scan;
                rx_apo=h.receive_apodization.data;
                rx_propagation_distance=h.receive_apodization.propagation_distance;
                
                % precalculate transmit apodization according to 10.1109/TUFFC.2015.007183
                % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
                h.transmit_apodization.sequence=h.channel_data.sequence(n_wave);
                h.transmit_apodization.scan=current_scan;
                tx_apo=h.transmit_apodization.data;
                
                % create an intermediate beamformed data class
                inter_dataset(n_wave)=uff.beamformed_data();
                inter_dataset(n_wave).scan=current_scan;
                inter_dataset(n_wave).wave=h.channel_data.sequence(n_wave);
                inter_dataset(n_wave).data=zeros(current_scan.N_pixels,h.channel_data.N_elements,h.channel_data.N_frames);
                temp_delayed = zeros(current_scan.N_pixels,h.channel_data.N_elements); %Using a temporary buffer is much faster
                
                % transmit delay
                if ~isinf(h.channel_data.sequence(n_wave).source.distance)
                    % point sources
                    TF=(-1).^(current_scan.z<h.channel_data.sequence(n_wave).source.z).*sqrt((h.channel_data.sequence(n_wave).source.x-current_scan.x).^2+(h.channel_data.sequence(n_wave).source.y-current_scan.y).^2+(h.channel_data.sequence(n_wave).source.z-current_scan.z).^2);
                    % add distance from source to origin
                    TF=TF+sign(cos(h.channel_data.sequence(n_wave).source.azimuth)).*h.channel_data.sequence(n_wave).source.distance;
                else
                    % plane waves
                    TF=current_scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+current_scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
                end
                
                % frame loop
                for n_frame=1:h.channel_data.N_frames
 
                    % receive loop
                    for nrx=1:h.channel_data.N_elements
                        
                        % Update waitbar
                        tools.workbar(counter/total_count,sprintf('%s (%s)',name,version),'USTB');
                        counter = counter+1;
                        
                        % receive delay
                        RF=sqrt((h.channel_data.probe.x(nrx)-current_scan.x).^2+(h.channel_data.probe.y(nrx)-current_scan.y).^2+(h.channel_data.probe.z(nrx)-current_scan.z).^2);
                        
                        % total delay
                        delay=(RF+TF)/h.channel_data.sound_speed;
                        
                        % check whether is IQ or RF data
                        if(w0>eps)
                            % phase correction factor
                            phase_shift=exp(1i.*w0*delay);
                            % data
                            data=h.channel_data.data(:,nrx,n_wave,n_frame);
                        else
                            % phase correction factor
                            phase_shift=1;
                            % data
                            data=hilbert(h.channel_data.data(:,nrx,n_wave,n_frame));
                        end
                        
                        %Using a temporary buffer is much faster
                        temp_delayed(:,nrx) = tx_apo.*rx_apo(:,nrx).*phase_shift.*interp1(h.channel_data.time,data,delay,'linear',0);
                    end
                    
                    inter_dataset(n_wave).data(:,:,n_frame) = temp_delayed;    
                    postprocess.bmf = h;
                    inter_dataset(n_wave) = postprocess.go(inter_dataset(n_wave));                
                    %inter_dataset(n_wave).data=bsxfun(@times,inter_dataset(n_wave).data,exp(-j*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
                end
            end
        end
