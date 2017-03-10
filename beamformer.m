classdef beamformer
%BEAMFORMER   beamformer definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/28 $

    %% public properties
    properties  (SetAccess = public)
        raw_data             % RAW_DATA class
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
            
            h.receive_apodization=apodization();  % APODIZATION class
            h.transmit_apodization=apodization(); % APODIZATION class
        end
    end
    
    %% set methods
    methods  
        function out_dataset=go(h,postprocess)
            
            % checking we have all we need
            assert(numel(h.raw_data)>0,'The RAW_DATA parameter is not set.');
            assert(numel(h.scan)>0,'The SCAN parameter is not set.');

            % modulation frequency
            w0=2*pi*h.raw_data.modulation_frequency;

            %% beamforming 
            wb=waitbar(0,sprintf('USTB General Beamformer (%s)',h.version));
            set(wb,'Name','USTB');
            for n_wave=1:numel(h.raw_data.sequence)
                waitbar(n_wave/numel(h.raw_data.sequence));
                
                % support multiple or single scans with the same code
                if numel(h.scan)==1 
                    current_scan=h.scan; 
                else
                    current_scan=h.scan(n_wave); 
                end

                % precalculate receive apodization
                h.receive_apodization.probe=h.raw_data.probe;
                h.receive_apodization.scan=current_scan;
                rx_apo=h.receive_apodization.data;
                rx_propagation_distance=h.receive_apodization.propagation_distance;
                
                % precalculate transmit apodization according to 10.1109/TUFFC.2015.007183 
                % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
                %h.transmit_apodization.probe=probe();
                h.transmit_apodization.sequence=h.raw_data.sequence(n_wave);
                h.transmit_apodization.scan=current_scan;
                tx_apo=h.transmit_apodization.data;
                
                % create an intermediate beamformed data class
                inter_dataset(n_wave)=beamformed_data();
                inter_dataset(n_wave).scan=current_scan; 
                inter_dataset(n_wave).wave=h.raw_data.sequence(n_wave);
                inter_dataset(n_wave).data=zeros(current_scan.N_pixels,1);
                
                % transmit delay
                if ~isinf(h.raw_data.sequence(n_wave).source.distance)
                    % point sources
                    TF=sqrt((h.raw_data.sequence(n_wave).source.x-current_scan.x).^2+(h.raw_data.sequence(n_wave).source.y-current_scan.y).^2+(h.raw_data.sequence(n_wave).source.z-current_scan.z).^2);
                    % add distance from source to origin
                    TF=TF+sign(cos(h.raw_data.sequence(n_wave).source.azimuth)).*h.raw_data.sequence(n_wave).source.distance;
                else
                    % plane waves
                    TF=h.scan.z*cos(h.raw_data.sequence(n_wave).source.azimuth)*cos(h.raw_data.sequence(n_wave).source.elevation)+h.scan.x*sin(h.raw_data.sequence(n_wave).source.azimuth)*cos(h.raw_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.raw_data.sequence(n_wave).source.elevation);
                end

                % receive loop
                for nrx=1:h.raw_data.probe.N_elements
                    
                    % receive delay
                    RF=sqrt((h.raw_data.probe.x(nrx)-current_scan.x).^2+(h.raw_data.probe.y(nrx)-current_scan.y).^2+(h.raw_data.probe.z(nrx)-current_scan.z).^2);

                    % total delay
                    delay=(RF+TF)/h.raw_data.sound_speed;

                    % phase correction factor (IQ)
                    phase_shift=exp(1i.*w0*delay);
               
                    % beamformed signal
                    inter_dataset(n_wave).data=inter_dataset(n_wave).data+tx_apo.*rx_apo(:,nrx).*phase_shift.*interp1(h.raw_data.time,h.raw_data.data(:,nrx,n_wave),delay,'linear',0);
                end
                
                % assign phase according to 2 times the receive propagation distance
                inter_dataset(n_wave).data=bsxfun(@times,inter_dataset(n_wave).data,exp(-j*w0*2*rx_propagation_distance/h.raw_data.sound_speed));
            end
            close(wb);

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
            assert(strcmp(class(in_pulse),'pulse'), 'The input is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.receive_apodization(h,in_apodization)
            assert(strcmp(class(in_apodization),'apodization'), 'The input is not a APODIZATION class. Check HELP APODIZATION.');
            h.receive_apodization=in_apodization;
        end        
        function h=set.transmit_apodization(h,in_apodization)
            assert(strcmp(class(in_apodization),'apodization'), 'The input is not a APODIZATION class. Check HELP APODIZATION.');
            h.transmit_apodization=in_apodization;
        end           
    end
    
    
end