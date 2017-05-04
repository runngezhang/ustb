classdef das_cc_standalone < beamformer
%DAS_CC_STANDALONE   Example of DAS and coherent compounding beamformer 
% using standalone MATLAB code
%
%   See also BEAMFORMER

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Date: 2017/03/10$

    %% constructor
    methods (Access = public)
        function h=das_cc_standalone(bmf)
            %DAS_CC_STANDALONE   Constructor of das_cc_standalone class
            %
            %   Syntax:
            %   h = das_cc_standalone()
            %
            %   See also BEAMFORMER                     
           
            if nargin<1
                bmf=[]; 
            end
            
            h=h@beamformer(bmf);
            
            % Logistics
            h.name='DAS Cohe Comp standalone';   
            h.reference= 'www.ustb.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};    
            h.version='v1.0.1';     
        end
    end
    
    
    %% go method
    methods 
        function out_data = go(h)
        
            % modulation frequency
            w0=2*pi*h.channel_data.modulation_frequency;

            % create output beamformed_data class
            out_data=uff.beamformed_data();
            out_data.scan=h.scan;
            out_data.data=zeros(h.scan.N_pixels,h.channel_data.N_frames);
            
            % wave loop
            tools.workbar();
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

                % receive loop
                for nrx=1:h.channel_data.N_elements
                    %waitbar(((n_wave-1)*h.channel_data.N_elements+nrx)/numel(h.channel_data.sequence)/h.channel_data.N_elements,wb);
                    tools.workbar(((n_wave-1)*h.channel_data.N_elements+nrx)/numel(h.channel_data.sequence)/h.channel_data.N_elements,sprintf('%s (%s)',h.name,h.version),'USTB');
                    % receive delay
                    RF=sqrt((h.channel_data.probe.x(nrx)-current_scan.x).^2+(h.channel_data.probe.y(nrx)-current_scan.y).^2+(h.channel_data.probe.z(nrx)-current_scan.z).^2);

                    % total delay
                    delay=(RF+TF)/h.channel_data.sound_speed;

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
                        out_data.data(:,n_frame)=out_data.data(:,n_frame)+tx_apo.*rx_apo(:,nrx).*phase_shift.*interp1(h.channel_data.time,data,delay,'linear',0);
                    end
                end

                % assign phase according to 2 times the receive propagation distance
                out_data.data=bsxfun(@times,out_data.data,exp(-j*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
            end
            tools.workbar(1);    
        end
    end    
end