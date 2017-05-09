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
        N_active_elements          % number of active transducers on receive
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
    
    % UFF variables
    properties (SetAccess = private)
        uff_version = 'v1.0.0';
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
                    eval(sprintf('h.%s = object.%s;',property_name,property_name));
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
    
    %% Ultrasound File Format
    methods (Access = public)    
        function write(h,filename)
            
            %-- Writes all the information to Ultrasound File Format (UFF)
            %-- Syntax:
            %-- write_file_hdf5(file_name)
            %-- file_name: Name of the hdf5 file
            
            
            open_write(filename);                                                % open UFF for writting
            location=uff.create_group(filename,[],'channel_data',h.uff_version); % create group
            
            % dump properties
            uff.append(filename, location, 'sound_speed', h.sound_speed);
            uff.append(filename, location, 'initial_time', h.initial_time);
            uff.append(filename, location, 'sampling_frequency', h.sampling_frequency);
            uff.append(filename, location, 'modulation_frequency', h.modulation_frequency);
            uff.append(filename, location, 'data', h.data);
            %uff.append(filename, location, 'probe', h.probe);

        end
       

%         function read_file(h,filename)
%             
%             %-- Reads all the information from a mat or hdf5 file
%             %-- Syntax:
%             %-- read_file(file_name)
%             %-- file_name: Name of the mat or hdf5 file
%             
%             [pathstr, name, ext] = fileparts(filename); 
%             switch ext
%                 case '.mat'
%                     h.read_file_mat(filename);
%                 case '.hdf5'
%                     h.read_file_hdf5(filename);
%                 otherwise
%                     error('Unknown signal format!');
%             end
%             
%         end
       
        
%         function write_file(h,filename)
%             
%             %-- Write all the information into a mat or hdf5 file
%             %-- Syntax:
%             %-- write_file(file_name)
%             %-- file_name: Name of the mat or hdf5 file
%             
%             [pathstr, name, ext] = fileparts(filename); 
%             switch ext
%                 case '.mat'
%                     h.write_file_mat(filename);
%                 case '.hdf5'
%                     h.write_file_hdf5(filename);
%                 otherwise
%                     error('Unknown signal format!');
%             end
%             
%         end
                
        
%         function read_file_hdf5(h,filename)
% 
%             %-- Reads all the information from a HUFF (HDF5 Ultrasound File Format) file
%             %-- Syntax:
%             %-- read_file_hdf5(file_name)
%             %-- file_name: Name of the hdf5 file
%             
%             %-- read US metagroup
%             info = h5info(filename,'/US');
% 
%             %-- read the groups in the metagroup
%             for n=1:length(info.Groups)
%                 location=info.Groups(n).Name;
%                 dstype=h5readatt(filename,location,'type');
%                 if strcmp(dstype,'US')
%                     subtype=h5readatt(filename,location,'subtype');
%                     if strcmp(subtype{1},'CPW')
% 
%                         %-- subtype
%                         dataset_subtype=h5readatt(filename,location,'subtype');
%                         assert(strcmp(dataset_subtype,'CPW'),'Only CPWC us_dataset are supported!');
%                         
%                         %-- read signal format 
%                         signal_format=h5readatt(filename,location,'signal_format');
%                         
%                         %-- read modulation frequency
%                         h.modulation_frequency=h5read(filename,[location '/modulation_frequency']);
% 
%                         %-- check format
%                         switch(signal_format{1})
%                             case 'RF'
%                                 assert(h.modulation_frequency==0,'RF dataset cannot have a modulation frequency');
%                             case 'IQ'
%                                 assert(h.modulation_frequency>0,'IQ dataset cannot have a null modulation frequency');
%                             otherwise
%                                 error('Unknown signal format!');
%                         end
% 
%                         %-- Attributes
%                         %-- read name
%                         a = h5readatt(filename,location,'name'); h.name=a{1}(1:end-1);
% 
%                         %-- read date
%                         a = h5readatt(filename,location,'creation_date'); h.creation_date=a{1}(1:end-1);
%                         
%                         %-- read version
%                         a = h5readatt(filename,location,'version'); h.version=a{1}(1:end-1);
% 
%                         %-- Data
%                         %-- read speed of sound
%                         h.c0 = h5read(filename,[location '/sound_speed']);
% 
%                         %-- read initial_time
%                         h.initial_time = h5read(filename,[location '/initial_time']);
% 
%                         %-- read sampling_frequency
%                         h.sampling_frequency = h5read(filename,[location '/sampling_frequency']);
% 
%                         %-- read sampling_frequency
%                         h.PRF = h5read(filename,[location '/PRF']);
% 
%                         %-- read transducer geometry
%                         h.probe_geometry = h5read(filename,[location '/probe_geometry']);
% 
%                         %-- read angles
%                         h.angles = h5read(filename,[location '/angles']);
%                         
%                         %-- read data
%                         real_part = h5read(filename,[location '/data/real']);
%                         imag_part = h5read(filename,[location '/data/imag']);
%                         h.data = real_part+1i*imag_part;
%                     end
%                 end
%             end
%         end

    end
end