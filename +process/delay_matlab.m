classdef delay_matlab < process
%DELAY_MATLAB   Matlab implementation of the Delay step of general beamformed
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/05/04$

    %% constructor
    methods (Access = public)
        function h=delay_matlab()
            h.name='Delay USTB General Beamformer MATLAB';   
            h.reference= 'www.ustb.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};    
            h.version='v1.0.3';
        end
    end

    methods
        function out_data=go(h)

            % modulation frequency
            w0=2*pi*h.channel_data.modulation_frequency;

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
                    TF=TF+sign(cos(h.channel_data.sequence(n_wave).source.azimuth)).*h.channel_data.sequence(n_wave).source.distance;
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
                for nrx=1:h.channel_data.N_elements
                    % progress bar
                    n=(n_wave-1)*h.channel_data.N_elements+nrx;
                    if mod(n,round(N/100))==1
                        tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                    end
                   
                    % create beamformed data class
                    out_data(nrx,n_wave)=uff.beamformed_data();
                    out_data(nrx,n_wave).scan=current_scan;
                    out_data(nrx,n_wave).wave=h.channel_data.sequence(n_wave);
                    out_data(nrx,n_wave).data=zeros(current_scan.N_pixels,h.channel_data.N_frames);

                    % total delay
                    delay=(RF(:,nrx)+TF)/h.channel_data.sound_speed;

                    for n_frame=1:h.channel_data.N_frames

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

                        % beamformed signal
                        out_data(nrx,n_wave).data(:,n_frame)=tx_apo.*rx_apo(:,nrx).*phase_shift.*interp1(h.channel_data.time,data,delay,'linear',0);
                    end
                    
                    % assign phase according to 2 times the receive propagation distance
                    out_data(nrx,n_wave).data=bsxfun(@times,out_data(nrx,n_wave).data,exp(-j*w0*2*rx_propagation_distance/h.channel_data.sound_speed));                    
                end
            end
            tools.workbar(1);
        end
    end
end
