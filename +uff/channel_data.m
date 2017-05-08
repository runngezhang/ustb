classdef channel_data < handle
    %CHANNEL_DATA   CHANNEL_DATA definition. Children of HANDLE class
    %
    %   See also CHANNEL_DATA/CHANNEL_DATA, BEAMFORMED_DATA, PHANTOM, PROBE
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %   $Date: 2017/02/24 $
    
    %% compulsory properties
    properties  (SetAccess = public)
        sampling_frequency         % sampling frequency [Hz]
        initial_time               % time of the initial sample [s]
        sound_speed     =1540      % reference sound speed [m/s]
        sequence                   % collection of WAVE classes
        probe                      % PROBE class
        data                       % data
        modulation_frequency       % modulation frequency [Hz]
    end
    
    %% optional properties
    properties  (SetAccess = public)
        phantom                    % PHANTOM class [optional]
        pulse                      % PULSE class [optional]
        PRF                        % pulse repetition frequency [Hz]
    end
    
    %% dependent properties
    properties  (Dependent)
        N_samples          % number of samples in the data
        N_elements         % number of elements in the probe
        N_channels         % number of elements in the probe
        N_waves            % number of transmitted waves
        N_frames           % number of frames
        time
        lambda             % wavelength [m]
    end
    
    
    %% constructor
    methods (Access = public)
        function h=channel_data()
            %CHANNEL_DATA   Constructor of CHANNEL_DATA class
            %
            %   Syntax:
            %   h = channel_data()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE
            
            h.modulation_frequency=0;
            h.sampling_frequency=0;
        end
    end
    
    %% copy
    methods (Access = public)
        function copy(h,object)
            %COPY    Copy the values from another channel_data
            %
            %   Syntax:
            %   COPY(object)
            %       object       Instance of a channel_data class
            %
            %   See also SCAN, WAVE, SOURCE
            assert(isa(object,class(h)),'Class of the input object is not identical');
            
            % we copy all non-dependent public properties
            list_properties=properties(object);
            for n=1:numel(list_properties)
                property_name=list_properties{n};
                mp = findprop(h,property_name);
                if strcmp(mp.GetAccess,'public')&&~mp.Dependent
                    eval(sprintf('h.%s = object.%s',property_name,property_name));
                end
            end
        end
    end
    
    %% plot methods
    methods
        function plot(h,n_beam)
            if nargin<2
                n_beam=1;
            end
            
            figure;
            if abs(h.modulation_frequency)>eps
                subplot(1,2,1);
                imagesc(1:h.N_elements,h.time*1e6,real(h.data(:,:,n_beam))); grid on; axis tight;
                xlabel('Channel');
                ylabel('time [\mus]');
                set(gca,'fontsize',14);
                title(sprintf('Real Part - Beam %d',n_beam));
                subplot(1,2,2);                
                imagesc(1:h.N_elements,h.time*1e6,imag(h.data(:,:,n_beam))); grid on; axis tight;
                xlabel('Channel');
                ylabel('time [\mus]');
                set(gca,'fontsize',14);
                title(sprintf('Imaginary Part - Beam %d',n_beam));
            else
                imagesc(1:h.N_elements,h.time*1e6,h.data(:,:,n_beam)); grid on; axis tight;
                xlabel('Channel');
                ylabel('time [\mus]');
                set(gca,'fontsize',14);
                title(sprintf('Beam %d',n_beam));
            end
        end
    end
    
    %% set methods
    methods
        function h=set.phantom(h,in_phantom)
            assert(isa(in_phantom,'uff.phantom'), 'The _phantom_ is not a PHANTOM class. Check HELP PHANTOM.');
            h.phantom=in_phantom;
        end
        function h=set.pulse(h,in_pulse)
            assert(isa(in_pulse,'uff.pulse'), 'The pulse is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.probe(h,in_probe)
            assert(isa(in_probe,'uff.probe'), 'The probe is not a PROBE class. Check HELP PROBE.');
            h.probe=in_probe;
        end
        function h=set.sequence(h,in_sequence)
            assert(isa(in_sequence,'uff.wave'), 'The sequence is not a WAVE class. Check HELP WAVE.');
            h.sequence=in_sequence;
        end
        function h=set.sampling_frequency(h,in_sampling_frequency)
            assert(numel(in_sampling_frequency)==1, 'The sampling frequency must be a scalar');
            h.sampling_frequency=in_sampling_frequency;
        end
        function h=set.modulation_frequency(h,in_modulation_frequency)
            assert(numel(in_modulation_frequency)==1, 'The sampling frequency must be a scalar');
            h.modulation_frequency=in_modulation_frequency;
        end
        function h=set.sound_speed(h,in_sound_speed)
            assert(numel(in_sound_speed)==1, 'The sound speed must be a scalar');
            h.sound_speed=in_sound_speed;
        end
        function h=set.initial_time(h,in_initial_time)
            assert(numel(in_initial_time)==1, 'The initial time must be a scalar');
            h.initial_time=in_initial_time;
        end
        function h=set.data(h,in_data)
            % checking needed inputs
            assert(~isempty(h.probe), 'The probe structure must be set before inserting the data.');
            assert(~isempty(h.sequence), 'The sequence structure must be set before inserting the data.');
            assert(~isempty(h.sampling_frequency), 'The sampling_frequency must be set before inserting the data.');
            assert(~isempty(h.initial_time), 'The initial_time must be set before inserting the data.');
            
            assert(size(in_data,2)==h.N_elements, 'The number of elements in the probe does not match the channels in the inserted data (2nd dimension).');
            assert(size(in_data,3)==h.N_waves, 'The number of waves in the sequence does not match the waves in the inserted data (3th dimension).');
            
            h.data=in_data;
        end
        function h=set.PRF(h,in_PRF)
            if ~isempty(in_PRF)
                assert(numel(in_PRF)==1, 'The PRF must be a scalar');
                h.PRF=in_PRF;
            end
        end
    end
    
    %% get methods
    methods
        function value=get.N_elements(h)
            value=h.probe.N_elements;
        end
        function value=get.N_channels(h)
            value=h.probe.N_elements;
        end
        function value=get.N_samples(h)
            value=size(h.data,1);
        end
        function value=get.N_waves(h)
            value=numel(h.sequence);
        end
        function value=get.N_frames(h)
            value=size(h.data,4);
        end
        function value=get.time(h)
            value=(h.initial_time+(0:h.N_samples-1)/h.sampling_frequency).';
        end
        function value=get.lambda(h)
            assert(~isempty(h.sound_speed),'You need to set the channel_data.sound_speed')
            assert(~isempty(h.pulse),'You need to set the pulse and the pulse center frequency.')
            value = h.sound_speed/h.pulse.center_frequency;
        end
    end
end