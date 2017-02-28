classdef simulator
%simulator   Simulator definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/24 $

    %% public properties
    properties  (SetAccess = public)
        phantom             % phantom class
        pulse               % pulse class
        probe               % probe class
        sequence            % collection of wave classes
        sampling_frequency  % sampling frequency [Hz]
    end
    
    %% dependent properties
    properties  (Dependent)   
        N_elements         % number of elements in the probe
        N_points           % number of points in the phantom
        N_beams            % number of beams 
    end
    
    %% private properties
    properties  (Access = private)   
        version='v1.0.2';  % simulator version
    end
    
    %% constructor
    methods (Access = public)
        function h=simulator()
            %simulator   Constructor of simulator class
            %
            %   Syntax:
            %   h = simulator()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE                      
            
        end
    end
    
    %% set methods
    methods  
        function out_dataset=go(h)
            disp(sprintf('USTB Fresnel linear impulse response simulator (%s)',h.version));
            disp('---------------------------------------------------------------');
            
            %% checking we have all we need
            assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(numel(h.phantom)>0,'The PHANTOM parameter is not set.');
            assert(numel(h.pulse)>0,'The PULSE parameter is not set.');
            assert(numel(h.sequence)>0,'The SEQUENCE parameter is not set.');
            assert(numel(h.sampling_frequency)>0,'The SAMPLING_FREQUENCY parameter is not set.');
            
            % checking number of elements
            assert(h.probe.N_elements==h.sequence(1).N_elements,'Mismatch in the number of elements in probe and the size of delay and apodization vectors in beam');
            
            %% time vector
            max_range=0;
            min_range=Inf;
            for n_p=1:h.N_points
                max_range=max([max_range; sqrt(sum((ones(h.N_elements,1)*h.phantom.points(n_p,1:3)-h.probe.geometry(:,1:3)).^2,2))]);
                min_range=min([min_range; sqrt(sum((ones(h.N_elements,1)*h.phantom.points(n_p,1:3)-h.probe.geometry(:,1:3)).^2,2))]);
            end
            time_1w=(min_range/h.phantom.sound_speed-8/h.pulse.center_frequency/h.pulse.fractional_bandwidth):(1/h.sampling_frequency):(max_range/h.phantom.sound_speed+4/h.pulse.center_frequency/h.pulse.fractional_bandwidth);                                                  % time vector [s]
            time_2w=(2*min_range/h.phantom.sound_speed-8/h.pulse.center_frequency/h.pulse.fractional_bandwidth):(1/h.sampling_frequency):(2*max_range/h.phantom.sound_speed+4/h.pulse.center_frequency/h.pulse.fractional_bandwidth);                                                  % time vector [s]
            N_samples=numel(time_2w);                                                                              % number of time samples

            %% wavenumber
            k=2*pi*h.pulse.center_frequency/h.phantom.sound_speed;
   
            %% unwrapping the signal
            focusing_delay=zeros(h.N_elements,1,h.N_beams);
            apodization=zeros(h.N_elements,1,h.N_beams);
            for n_beam=1:h.N_beams 
                focusing_delay(:,1,n_beam)=h.sequence(n_beam).delay;
                apodization(:,1,n_beam)=h.sequence(n_beam).apodization;
            end
                  
            %% points loop
            data=zeros(N_samples,h.N_elements,h.N_beams);
            wb = waitbar(0,'Simulating');
            for n_p=1:h.N_points
                % waitbar
                waitbar(n_p/h.N_points,wb,sprintf('Simulating point %d of %d',n_p,h.N_points));
                
                % computing geometry relations to the point
                theta     = atan2(h.phantom.x(n_p)-h.probe.x, h.phantom.z(n_p)-h.probe.z)-h.probe.theta;
                phi       = atan2(h.phantom.y(n_p)-h.probe.y, h.phantom.z(n_p)-h.probe.z)-h.probe.phi;
                distance  = sqrt(sum((h.probe.geometry(:,1:3)-ones(h.N_elements,1)*h.phantom.points(n_p,1:3)).^2,2));

                % directivity between probe and the point
                directivity = sinc(k*h.probe.width/2.*sin(theta).*cos(phi)/pi).*sinc(k*h.probe.height/2.*sin(theta).*sin(phi)/pi);
                % delay between probe and the point
                propagation_delay = distance/h.phantom.sound_speed;
                    
%                 % the beam loop
%                 for n_beam=1:h.N_beams 
%                     
%                     % computing the transmit signal 
%                     delayed_time=ones(h.N_elements,1)*time-(propagation_delay+h.sequence(n_beam).delay)*ones(1,N_samples);                
%                     transmit_signal=sum(bsxfun(@times,h.pulse.signal(delayed_time),h.sequence(n_beam).apodization.*directivity./(4*pi*distance)),1);  
% 
%                     % computing the receive signal
%                     delayed_time=ones(h.N_elements,1)*time-propagation_delay*ones(1,N_samples);                
%                     data(:,:,n_beam)=data(:,:,n_beam)+bsxfun(@times,interp1(time,transmit_signal,delayed_time,'linear',0),10.^(-h.phantom.alpha*(distance/1e-2)*h.pulse.center_frequency).*directivity./(4*pi*distance)).';                      
%                 end
                
                %% all in one
                delayed_time_1w=ones(h.N_elements,1)*time_1w-(propagation_delay)*ones(size(time_1w));  
                delayed_time_2w=ones(h.N_elements,1)*time_2w-(propagation_delay)*ones(1,N_samples);  
                DDT=repmat(delayed_time_1w,[1 1 h.N_beams]);
                mf=h.phantom.Gamma(n_p).*10.^(-h.phantom.alpha*(distance/1e-2)*(h.pulse.center_frequency/1e6)/20).*directivity./(4*pi*distance);
                AMF=bsxfun(@times,apodization,mf);
                transmit_signal=squeeze(sum(bsxfun(@times,h.pulse.signal(bsxfun(@minus,DDT,focusing_delay)),AMF),1));
                data=data+permute(bsxfun(@times,interp1(time_1w,transmit_signal,delayed_time_2w,'linear',0),mf),[2 1 3]);
            end
            delete(wb);
            
            % save the data into a RAW_DATA structure
            out_dataset=raw_data();
            out_dataset.probe=h.probe();
            out_dataset.pulse=h.pulse();
            out_dataset.phantom=h.phantom();
            out_dataset.sequence=h.sequence();
            out_dataset.sampling_frequency=h.sampling_frequency();
            out_dataset.sound_speed=h.phantom.sound_speed;
            out_dataset.initial_time=time_2w(1);
            out_dataset.data=data;
            
        end
    end
    
    %% set methods
    methods  
        function h=set.phantom(h,in_phantom)
            assert(strcmp(class(in_phantom),'phantom'), 'The phantom is not a PHANTOM class. Check HELP PHANTOM.');
            h.phantom=in_phantom;
        end
        function h=set.pulse(h,in_pulse)
            assert(strcmp(class(in_pulse),'pulse'), 'The pulse is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.probe(h,in_probe)
            assert(strcmp(class(in_probe),'probe'), 'The probe is not a PROBE class. Check HELP PROBE.');
            h.probe=in_probe;
        end
        function h=set.sequence(h,in_sequence)
            assert(strcmp(class(in_sequence),'wave'), 'The sequence members are not a WAVE class. Check HELP WAVE.');
            h.sequence=in_sequence;
        end
        function h=set.sampling_frequency(h,in_sampling_frequency)
            assert(numel(in_sampling_frequency)==1, 'The sampling frequency must be a escalar');
            h.sampling_frequency=in_sampling_frequency;
        end       
    end
    
    %% get methods
    methods  
        function value=get.N_elements(h)
            value=h.probe.N_elements;
        end
        function value=get.N_points(h)
            value=h.phantom.N_points;
        end
        function value=get.N_beams(h)
            value=numel(h.sequence);
        end
    end
    
end