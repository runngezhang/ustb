classdef das_mex < process
%DAS_MEX   Mex implementation of the Delay-and-Sum general beamformer
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%
%   $Last updated: 2017/05/04$

    %% constructor
    methods (Access = public)
        function h=das_mex()
            h.name='USTB General DAS Beamformer MEX';   
            h.reference= 'www.ustb.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};    
            h.version='v1.0.3';          
        end
    end
    
    %% go method
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

            % precalculating hilbert (if needed)
            data=h.channel_data.data;
            if ~(w0>eps)
                data=single(reshape(hilbert(h.channel_data.data(:,:)),size(h.channel_data.data)));
            end
            
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

                % create an intermediate beamformed data class
                out_data(1,n_wave)=uff.beamformed_data();
                out_data(1,n_wave).scan=current_scan;
                out_data(1,n_wave).wave=h.channel_data.sequence(n_wave);
                out_data(1,n_wave).data=zeros(current_scan.N_pixels,h.channel_data.N_frames);

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

                % total delay
                delay=bsxfun(@plus,RF,TF)./h.channel_data.sound_speed;
                
                % factor
                apodization_matrix=bsxfun(@times,tx_apo,rx_apo);%.*phase_shift;
                
                % das
                out_data(1,n_wave).data=mex.das_c(single(squeeze(data(:,:,n_wave,:))),single(h.channel_data.sampling_frequency),single(h.channel_data.initial_time),single(delay),single(apodization_matrix),single(h.channel_data.modulation_frequency));                

                % assign phase according to 2 times the receive propagation distance
                if(w0>eps)
                    out_data(1,n_wave).data=bsxfun(@times,out_data(1,n_wave).data,exp(-1i*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
                end
            end
            tools.workbar(1);
        end
    end
end
