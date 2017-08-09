classdef das_matlab < process
%DAS_MATLAB   Matlab implementation of the Delay-and-Sum general beamformed
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/07/09$

    %% constructor
    methods (Access = public)
        function h=das_matlab()
            h.name='USTB General DAS Beamformer MATLAB';   
            h.reference= 'www.ustb.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};    
            h.version='v1.0.7';          
        end
    end
    
    %% go method
    methods
        function out_data=go(h)

            % modulation frequency
            w0=2*pi*h.channel_data.modulation_frequency;

            % precalculate receive apodization
            h.receive_apodization.probe=h.channel_data.probe;
            h.receive_apodization.focus=h.scan(1);
            rx_apo=h.receive_apodization.data;
            rx_propagation_distance=h.receive_apodization.propagation_distance;
            
            % precalculate receive delay
            xm=bsxfun(@minus,h.channel_data.probe.x.',h.scan(1).x);
            ym=bsxfun(@minus,h.channel_data.probe.y.',h.scan(1).y);
            zm=bsxfun(@minus,h.channel_data.probe.z.',h.scan(1).z);
            RF=sqrt(xm.^2+ym.^2+zm.^2);

            % precalculating hilbert (if needed)
            data=h.channel_data.data;
            if ~(w0>eps)
                data=reshape(hilbert(h.channel_data.data(:,:)),size(h.channel_data.data));
            end
            
            % create beamformed data class
            out_data=uff.beamformed_data();
            out_data.scan=h.scan;
            N_pixels = 0; for n=1:length(h.scan) N_pixels = max([N_pixels h.scan(n).N_pixels]); end
            % out_data.sequence=h.channel_data.sequence; % not included by default
                           
            % auxiliary data
            aux_data=zeros(N_pixels,1,numel(h.channel_data.sequence),h.channel_data.N_frames);
 
            % wave loop
            tools.workbar();
            N=numel(h.channel_data.sequence)*h.channel_data.N_elements;
            for n_wave=1:numel(h.channel_data.sequence)

                % support multiple or single scans with the same code
                if numel(h.scan)==1
                    current_scan=h.scan;
                else
                    current_scan=h.scan(n_wave);
                end

                % calculate receive apodization for multiple scan
                if numel(h.scan)>1
                    h.receive_apodization.focus=current_scan;
                    rx_apo=h.receive_apodization.data;
                    rx_propagation_distance=h.receive_apodization.propagation_distance;
                end 
                
                % precalculate transmit apodization according to 10.1109/TUFFC.2015.007183
                % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
                h.transmit_apodization.sequence=h.channel_data.sequence(n_wave);
                h.transmit_apodization.focus=current_scan;
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
                
                % receive loop
                for n_rx=1:h.channel_data.N_elements
                    % progress bar
                    n=(n_wave-1)*h.channel_data.N_elements+n_rx;
                    if mod(n,round(N/100))==1
                        tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                    end
                    
                    % total delay
                    delay=(RF(:,n_rx)+TF)/h.channel_data.sound_speed;

                    for n_frame=1:h.channel_data.N_frames

                        % phase correction factor
                        if(w0>eps)
                            phase_shift=exp(1i.*w0*delay);
                        else
                            phase_shift=1;
                        end

                        % beamformed signal
                        aux_data(:,1,n_wave,n_frame)=aux_data(:,1,n_wave,n_frame)+tx_apo.*rx_apo(:,n_rx).*phase_shift.*interp1(h.channel_data.time,data(:,n_rx,n_wave,n_frame),delay,'linear',0);
                    end
                end

                % assign phase according to 2 times the receive propagation distance
                aux_data(:,1,n_wave,:)=bsxfun(@times,aux_data(:,1,n_wave,:),exp(-j*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
            end
            out_data.data=aux_data;
            tools.workbar(1);
        end
    end
end
