classdef dataset
%dataset   Dataset definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2016/09/01 $

    %% public properties
    properties  (SetAccess = public)
        sampling_frequency         % sampling frequency [Hz]
        initial_time               % time of the initial sample [s]
        probe                      % probe class
        pulse                      % pulse class
        phantom                    % phantom class
        sequence                   % collection of beam class
        data                       % data
    end
    
    %% dependent properties
    properties  (Dependent)   
        N_elements         % number of elements in the probe
        N_samples          % number of samples in the data
        N_beams            % number of beams 
        time
    end
    
    
    %% constructor
    methods (Access = public)
        function h=dataset()
            %dataset   Constructor of dataset class
            %
            %   Syntax:
            %   h = dataset()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE
                      
        end
    end
    
    %% plot methods
    methods
        function plot(h,n_beam)
           if nargin<2
               n_beam=1;
            end
            
          t=h.initial_time+(0:h.N_samples-1)/h.sampling_frequency;
          
          figure; 
          imagesc(1:h.N_elements,t*1e6,h.data(:,:,n_beam)); grid on; axis tight;
          xlabel('Channel');
          ylabel('time [\mus]');
          set(gca,'fontsize',14);
          title(sprintf('Beam %d',n_beam));
        end
    end
    
    %% set methods
    methods  
        function h=set.phantom(h,in_phantom)
            assert(strcmp(class(in_phantom),'phantom'), 'The phantom_ is not a PHANTOM class. Check HELP PHANTOM.');
            h.phantom=in_phantom;
        end
        function h=set.pulse(h,in_pulse)
            assert(strcmp(class(in_pulse),'pulse'), 'The pulse_ is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.probe(h,in_probe)
            assert(strcmp(class(in_probe),'probe'), 'The probe_ is not a PROBE class. Check HELP PROBE.');
            h.probe=in_probe;
        end
        function h=set.sequence(h,in_sequence)
            assert(strcmp(class(in_sequence),'beam'), 'The sequence_ is not a BEAM class. Check HELP BEAM.');
            h.sequence=in_sequence;
        end
        function h=set.sampling_frequency(h,in_sampling_frequency)
            assert(numel(in_sampling_frequency)==1, 'The sampling frequency must be a escalar');
            h.sampling_frequency=in_sampling_frequency;
        end 
        function h=set.initial_time(h,in_initial_time)
            assert(numel(in_initial_time)==1, 'The sampling frequency must be a escalar');
            h.initial_time=in_initial_time;
        end 
        function h=set.data(h,in_data)
            % checking needed inputs
            assert(~isempty(h.probe), 'The probe structure must be set before inserting the data.');
            assert(~isempty(h.pulse), 'The pulse structure must be set before inserting the data.');
            assert(~isempty(h.sequence), 'The sequence structure must be set before inserting the data.');
            assert(~isempty(h.sampling_frequency), 'The sampling_frequency must be set before inserting the data.');
            assert(~isempty(h.initial_time), 'The initial_time must be set before inserting the data.');
            
            assert(size(in_data,2)==h.N_elements, 'The N_elements in the probe does not match the channels in the inserted data (2nd dimension).');
            assert(size(in_data,3)==h.N_beams, 'The N_beams in the sequence does not match the beams in the inserted data (3th dimension).');
            
            h.data=in_data;
        end 
    end
    
    %% get methods
    methods  
        function value=get.N_elements(h)
            value=h.probe.N_elements;
        end
        function value=get.N_samples(h)
            value=size(h.data,1);
        end
        function value=get.N_beams(h)
            value=numel(h.sequence);
        end
        function value=get.time(h)
            value=h.initial_time+(0:h.N_samples-1)/h.sampling_frequency;
        end
    end
end