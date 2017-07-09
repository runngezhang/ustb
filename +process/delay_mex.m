classdef delay_mex < process
    %DELAY_MEX   Mex implementation of the Delay step of USTB general beamformer
    %
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %
    %   $Last updated: 2017/07/09$
    
    %% constructor
    methods (Access = public)
        function h=delay_mex()
            h.name='USTB Delay General Beamformer MEX';
            h.reference= 'www.ustb.no';
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};
            h.version='v1.0.7';
        end
    end
    
    %% go method
    methods
        function out_data=go(h)
            
            % Check if mex is properly setup on this computer
            if numel(mex.getCompilerConfigurations()) == 0 %Then no c/c++ compiler is set up.
                warning('No (c/c++)-compiler set up with MATLAB. Unable to use MEX, will default to MATLAB implementation.');
                matlab_delay = process.das_matlab();
                matlab_delay.channel_data = h.channel_data;
                matlab_delay.scan = h.scan;
                matlab_delay.transmit_apodization = h.transmit_apodization;
                matlab_delay.receive_apodization = h.receive_apodization;
                out_data = matlab_delay.go();
            else % MEX is ok!
                % modulation frequency
                w0=2*pi*h.channel_data.modulation_frequency;
                
                % constants
                sampling_frequency=single(h.channel_data.sampling_frequency);
                initial_time=single(h.channel_data.initial_time);
                modulation_frequency=single(h.channel_data.modulation_frequency);
                
                % precalculate receive apodization
                h.receive_apodization.probe=h.channel_data.probe;
                h.receive_apodization.scan=h.scan(1);
                rx_apo=h.receive_apodization.data;
                rx_propagation_distance=h.receive_apodization.propagation_distance;
                
                % precalculate receive delay
                xm=bsxfun(@minus,h.channel_data.probe.x.',h.scan(1).x);
                ym=bsxfun(@minus,h.channel_data.probe.y.',h.scan(1).y);
                zm=bsxfun(@minus,h.channel_data.probe.z.',h.scan(1).z);
                RF=sqrt(xm.^2+ym.^2+zm.^2);
                
                % precalculating hilbert (if needed)
                data=single(h.channel_data.data);
                if ~(w0>eps)
                    data=single(reshape(hilbert(h.channel_data.data(:,:)),size(h.channel_data.data)));
                end
                
                % create beamformed data class
                out_data=uff.beamformed_data();
                out_data.scan=h.scan;
                N_pixels = 0; for n=1:length(h.scan) N_pixels = max([N_pixels h.scan(n).N_pixels]); end
                %out_data.sequence=h.channel_data.sequence;
                
                % auxiliary data
                aux_data=zeros(N_pixels,h.channel_data.N_channels,numel(h.channel_data.sequence),h.channel_data.N_frames);
                
                % wave loop
                tools.workbar();
                N=numel(h.channel_data.sequence);
                for n_wave=1:numel(h.channel_data.sequence)
                    % progress bar
                    tools.workbar(n_wave/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                    
                    % support multiple or single scans with the same code
                    if numel(h.scan)==1
                        current_scan=h.scan;
                    else
                        current_scan=h.scan(n_wave);
                    end
                    
                    % calculate receive apodization for multiple scan
                    if numel(h.scan)>1
                        h.receive_apodization.scan=current_scan;
                        rx_apo=h.receive_apodization.data;
                        rx_propagation_distance=h.receive_apodization.propagation_distance;
                    end
                    
                    % precalculate transmit apodization according to 10.1109/TUFFC.2015.007183
                    % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
                    h.transmit_apodization.sequence=h.channel_data.sequence(n_wave);
                    h.transmit_apodization.scan=current_scan;
                    tx_apo=h.transmit_apodization.data;
                    
                    % transmit delay
                    if ~isinf(h.channel_data.sequence(n_wave).source.distance)
                        % point sources
                        TF=(-1).^(current_scan.z<h.channel_data.sequence(n_wave).source.z).*sqrt((h.channel_data.sequence(n_wave).source.x-current_scan.x).^2+(h.channel_data.sequence(n_wave).source.y-current_scan.y).^2+(h.channel_data.sequence(n_wave).source.z-current_scan.z).^2);
                        % add distance from source to origin
                        %OLD VERSION with rounding problem: TF=TF+sign(cos(h.channel_data.sequence(n_wave).source.azimuth)).*h.channel_data.sequence(n_wave).source.distance;
                        if (h.channel_data.sequence(n_wave).source.z<-1e-3)
                            TF=TF-h.channel_data.sequence(n_wave).source.distance;
                        else
                            TF=TF+h.channel_data.sequence(n_wave).source.distance;                        
                        end
                    else
                        % plane waves
                        TF=current_scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+current_scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
                    end
                    
                    % calculate receive delay for multiple scan
                    if numel(h.scan)>1
                        xm=bsxfun(@minus,h.channel_data.probe.x.',current_scan.x);
                        ym=bsxfun(@minus,h.channel_data.probe.y.',current_scan.y);
                        zm=bsxfun(@minus,h.channel_data.probe.z.',current_scan.z);
                        RF=sqrt(xm.^2+ym.^2+zm.^2);
                    end
                    
                    % total delay
                    delay=single(bsxfun(@plus,RF,TF)./h.channel_data.sound_speed);
                    
                    % factor
                    apodization_matrix=single(bsxfun(@times,tx_apo,rx_apo));
                    
                    % delay
                    aux_data(:,:,n_wave,:)=mex.delay_c(data(:,:,n_wave,:),sampling_frequency,initial_time,delay,apodization_matrix,modulation_frequency);
                    
                    % assign phase according to 2 times the receive propagation distance
                    if(w0>eps)
                        aux_data(:,:,n_wave,:)=bsxfun(@times,aux_data(:,:,n_wave,:),exp(-1i*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
                    end
                end
                out_data.data = aux_data;
                tools.workbar(1);
            end
        end
    end
end
