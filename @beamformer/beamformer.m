classdef beamformer
%BEAMFORMER   beamformer definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%            Stine Myren Hverven <stinemhv@ifi.uio.no>
%
%   $Date: 2017/03/10$

    %% public properties
    properties  (SetAccess = public)
        channel_data         % CHANNEL_DATA class
        receive_apodization  % APODIZATION class
        transmit_apodization % APODIZATION class
        scan                 % collection of SCAN classes
    end
    
    %% optional properties
    properties  (SetAccess = public)
        pulse                % PULSE class
    end
    
    %% private properties
    properties  (Access = private)   
        version='v1.0.3';  % beamformer version
    end
    
    %% constructor
    methods (Access = public)
        function h=beamformer()
            %beamformer   Constructor of beamformer class
            %
            %   Syntax:
            %   h = beamformer()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE                      
            
            h.receive_apodization=huff.apodization();  % APODIZATION class
            h.transmit_apodization=huff.apodization(); % APODIZATION class
        end
    end
    
    %% set methods
    methods  
        function out_dataset=go(h,postprocess)
            
            % checking we have all we need
            assert(numel(h.channel_data)>0,'The channel_data parameter is not set.');
            assert(numel(h.scan)>0,'The SCAN parameter is not set.');

            % modulation frequency
            w0=2*pi*h.channel_data.modulation_frequency;

            %% beamforming 
            %wb=waitbar(0,sprintf('USTB General Beamformer (%s)',h.version));
            %set(wb,'Name','USTB');
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
                inter_dataset(n_wave)=huff.beamformed_data();
                inter_dataset(n_wave).scan=current_scan; 
                inter_dataset(n_wave).wave=h.channel_data.sequence(n_wave);
                inter_dataset(n_wave).data=zeros(current_scan.N_pixels,1,h.channel_data.N_frames);

                % transmit delay
                if ~isinf(h.channel_data.sequence(n_wave).source.distance)
                    % point sources
                    %TF=sqrt((h.channel_data.sequence(n_wave).source.x-current_scan.x).^2+(h.channel_data.sequence(n_wave).source.y-current_scan.y).^2+(h.channel_data.sequence(n_wave).source.z-current_scan.z).^2);                    
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
                    tools.workbar(((n_wave-1)*h.channel_data.N_elements+nrx)/numel(h.channel_data.sequence)/h.channel_data.N_elements,sprintf('USTB General Beamformer (%s)',h.version),'USTB');
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
                        inter_dataset(n_wave).data(:,1,n_frame)=inter_dataset(n_wave).data(:,1,n_frame)+tx_apo.*rx_apo(:,nrx).*phase_shift.*interp1(h.channel_data.time,data,delay,'linear',0);
                    end
                end

                % assign phase according to 2 times the receive propagation distance
                inter_dataset(n_wave).data=bsxfun(@times,inter_dataset(n_wave).data,exp(-j*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
            end
            %close(wb);

            %% postprocess
            if nargin==1
                out_dataset=inter_dataset;
            else
                out_dataset=postprocess(inter_dataset);
            end
            
        end
    end
    
    %% set methods
    methods  
        function h=set.pulse(h,in_pulse)
            assert(isa(in_pulse,'huff.pulse'), 'The input is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.receive_apodization(h,in_apodization)
            assert(isa(in_apodization,'huff.apodization'), 'The input is not a APODIZATION class. Check HELP APODIZATION.');
            h.receive_apodization=in_apodization;
        end        
        function h=set.transmit_apodization(h,in_apodization)
            assert(isa(in_apodization,'huff.apodization'), 'The input is not a APODIZATION class. Check HELP APODIZATION.');
            h.transmit_apodization=in_apodization;
        end           
        function h=set.channel_data(h,in_channel_data)
            assert(isa(in_channel_data,'huff.channel_data'), 'The input is not a CHANNEL_DATA class. Check HELP CHANNEL_DATA.');
            h.channel_data=in_channel_data;
        end   
    end
    
    
end