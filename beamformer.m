classdef beamformer
%BEAMFORMER   beamformer definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/28 $

    %% public properties
    properties  (SetAccess = public)
        raw_data             % RAW_DATA class
        probe                % PROBE class
        receive_apodization  % APODIZATION class
        transmit_apodization % APODIZATION class
        sequence             % collection of WAVE classes
        scan                 % collection of SCAN classes
    end
    
    %% optional properties
    properties  (SetAccess = public)
        pulse                % PULSE class
    end
    
    %% private properties
    properties  (Access = private)   
        version='v1.0.2';  % beamformer version
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
        function out_dataset=go(h)
            disp(sprintf('USTB General beamformer (%s)',h.version));
            disp('---------------------------------------------------------------');
            
            % checking we have all we need
            assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(numel(h.sequence)>0,'The SEQUENCE parameter is not set.');
            assert(numel(h.scan)>0,'The SCAN parameter is not set.');

            % modulation frequency
            w0=2*pi*h.raw_data.modulation_frequency;

            % setting probe into apodization
            h.receive_apodization.probe=h.probe;
            h.transmit_apodization.probe=h.probe;
            
            % beamforming loop
            wb=waitbar(0,'Beamforming');
            for n_wave=1:numel(h.sequence)
                waitbar(n_wave/numel(h.sequence));
                
                % support multiple or single scans with the same code
                if numel(h.scan)==1 
                    current_scan=h.scan; 
                else
                    current_scan=h.scan(n_wave); 
                end

                % precalculate apodizations
                h.receive_apodization.scan=current_scan;
                h.transmit_apodization.scan=current_scan;
                rx_apo=h.receive_apodization.data;
                tx_apo=h.transmit_apodization.data;
                
                % create beamformed data class
                out_dataset(n_wave)=beamformed_data();
                out_dataset(n_wave).scan=current_scan; 
                out_dataset(n_wave).wave=h.sequence(n_wave);
                out_dataset(n_wave).data=zeros(current_scan.N_pixels,1);
                
                % transmit delay
                TF=sqrt((h.sequence(n_wave).source.x-current_scan.x).^2+(h.sequence(n_wave).source.y-current_scan.y).^2+(h.sequence(n_wave).source.z-current_scan.z).^2);

                % add distance from source to origin
                if ~isinf(h.sequence(n_wave).source.distance)
                    TF=TF+sign(cos(h.sequence(n_wave).source.azimuth)).*h.sequence(n_wave).source.distance;
                end

                % receive loop
                for nrx=1:h.probe.N_elements
                    
                    % receive delay
                    RF=sqrt((h.probe.x(nrx)-current_scan.x).^2+(h.probe.y(nrx)-current_scan.y).^2+(h.probe.z(nrx)-current_scan.z).^2);

                    % total delay
                    delay=(RF+TF)/h.raw_data.sound_speed;

                    % phase correction factor (IQ)
                    phase_shift=-exp(1i.*w0*delay);

                    % beamformed signal
                    out_dataset(n_wave).data=out_dataset(n_wave).data+rx_apo(:,nrx).*phase_shift.*interp1(h.raw_data.time,h.raw_data.data(:,nrx,n_wave),delay,'linear',0);
                end
            end
            close(wb);

        end
    end
    
    %% set methods
    methods  
        function h=set.pulse(h,in_pulse)
            assert(strcmp(class(in_pulse),'pulse'), 'The input is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.probe(h,in_probe)
            assert(strcmp(class(in_probe),'probe'), 'The input is not a PROBE class. Check HELP PROBE.');
            h.probe=in_probe;
        end
        function h=set.sequence(h,in_sequence)
            assert(strcmp(class(in_sequence),'wave'), 'The input is not a WAVE class. Check HELP WAVE.');
            h.sequence=in_sequence;
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