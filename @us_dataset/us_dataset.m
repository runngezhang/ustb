classdef us_dataset < handle
%US_DATASET    Ultrasound dataset.
%
%   US_DATASET is a superclass containing the properties and methods that are
%   common to a number of ultrasound datasets: Synthetic transmitaperture (sta), 
%   Coherent plane wave (cpw), Virtual source (vs). 
%
%   See also STA, CPW, VS

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $

    properties  (SetAccess = public)
        name=''             % String containing the name of the dataset
        creation_date       % String containing the date the dataset class was created
        format = E.signal_format.IQ  % Signal_format RF or IQ (enumerations.signal_format.RF, default=enumerations.signal_format.IQ)
        geom                % matrix M x 3 containing probe geometry [x, y, z] (m)
        data                % Matrix containing the numerical data. For acquisition
                            % datasets the matrix dimensions are [time_samples, channels, firings, frames]
                            % For image reconstruction datasets the matrix
                            % dimensions are [pixels, frames]
        time                % vector containing fast time (s)
        c0                  % reference speed of sound (m/s)
        center_frequency    % transmitted center frequency (Hz)
        modulation_frequency % value conatining the modulation frequency (Hz), only required for IQ format
    end
    
    properties  (SetAccess = protected)   
        frames                  % number of frames in the dataset        
        channels                % number of channels in the transducer
        firings                 % number of firings in the sequence (i.e. number of plane waves in CPWI, virtual sources in VSI, or elements in STAI)
        initial_time            % initial time (s)
        sampling_frequency      % sampling frequency (Hz)
        transmit_apodization    % matrix containing the apodization used for transmiting
        receive_apodization     % matrix containing the apodization used for receiving
    end
    
    %% constructor
    methods (Access = public)
        function h = us_dataset(name)
            %US_DATASET    Constructor of the DATASET class.
            %
            %   Syntax:
            %   US_DATASET(name) 
            %       name                    Name of the dataset
            %
            %   See also STA, CPW, VS
 
            if exist('name') h.name=name; end
            h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
        end
    end
    
    %% Relaunching method
    methods (Access = public)
        function sig=launch_implementation(h,r,imp)
            % Launches appropriate implementation
            switch(imp)
                case E.implementation.simple_matlab
                    sig=h.ir_simple_matlab(r);
                case E.implementation.optimized_matlab
                    sig=h.ir_optimized_matlab(r);
                case E.implementation.mex
                    sig=h.ir_mex(r);
                case E.implementation.mex_gpu
                    sig=h.ir_mex_gpu(r);
                otherwise
                    error('Selected implementation is not supported');
            end
        end
    end
    
    %% HUFF
    methods (Access = public)
        function huff_write(h,filename,location)
            %HUFF_WRITE    Dumps all the information of the ultrasound
            %dataset to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_WRITE(file_name,location)
            %       file_name                    Name of the hdf5 file
            %       location                   Name of the group destination
            %
            %   See also US_DATASET
            
            % Dataset type
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'US', 0); 
                H5T.enum_insert (filetype, 'SR', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'type', filetype, space, 'H5P_DEFAULT');
            H5A.write (attr, filetype, uint32(0));  % <--- US
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
                        
            % Signal format 
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'RF', 0); 
                H5T.enum_insert (filetype, 'IQ', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'signal_format', filetype, space, 'H5P_DEFAULT');
            switch(h.format)
                case E.signal_format.RF
                    H5A.write (attr, filetype, uint32(0));
                case E.signal_format.IQ
                    H5A.write (attr, filetype, uint32(1));
            end
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);

            % add name
            attr = h.name;
            attr_details.Name = 'name';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');

            % add date
            attr = h.creation_date;
            attr_details.Name = 'creation_date';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');

            % add speed of sound
            attr = h.c0;
            attr_details.Name = 'sound_speed';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');

            % add initial_time
            attr = h.time(1);
            attr_details.Name = 'initial_time';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');

            % add sampling_frequency
            attr = 1./mean(diff(h.time));
            attr_details.Name = 'sampling_frequency';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');

            % add modulation frequency
            if h.format==E.signal_format.IQ
                attr = h.modulation_frequency;
                attr_details.Name = 'modulation_frequency';
                attr_details.AttachedTo = location;
                attr_details.AttachType = 'group';
                hdf5write(filename, attr_details, attr, 'WriteMode', 'append');
            end

            % add data
            dset_details.Location = [location '/data'];
            dset_details.Name = 'real';
            hdf5write(filename, dset_details, real(h.data), 'WriteMode', 'append');
            dset_details.Name = 'imag';
            hdf5write(filename, dset_details, imag(h.data), 'WriteMode', 'append');
            
            % add transducer geometry
            dset_details.Location = [location '/transducer'];
            dset_details.Name = 'geom';
            hdf5write(filename, dset_details, h.geom, 'WriteMode', 'append');
        end
        
        function huff_read(h,filename,location)
            %HUFF_READ    Dumps all the information of the ultrasound
            %dataset to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_READ(file_name,location)
            %       file_name                    Name of the hdf5 file
            %       location                   Name of the group destination
            %
            %   See also US_DATASET
            
            % type
            dataset_type=h5readatt(filename,location,'type');
            assert(strcmp(dataset_type,'US'),'Ultrasound signal type does not match!');
                        
            % read signal format 
            signal_format=h5readatt(filename,location,'signal_format');
            switch(signal_format{1})
                case 'RF'
                    h.format=E.signal_format.RF;
                case 'IQ'
                    h.format=E.signal_format.IQ;
                otherwise
                    error('Unknown signal format!');
            end

            % read name
            h.name=h5readatt(filename,location,'name');

            % read date
            h.creation_date=h5readatt(filename,location,'creation_date');

            % read speed of sound
            h.c0=h5readatt(filename,location,'sound_speed');

            % read initial_time
            h.initial_time=h5readatt(filename,location,'initial_time');

            % read sampling_frequency
            h.sampling_frequency=h5readatt(filename,location,'sampling_frequency');
            
            % read modulation frequency
            if h.format==E.signal_format.IQ
                h.modulation_frequency=h5readatt(filename,location,'modulation_frequency');
            else
                h.modulation_frequency=0;
            end

            % read data
            real_part=h5read(filename,[location '/data/real']);
            imag_part=h5read(filename,[location '/data/imag']);
            h.data=real_part+1i*imag_part;
            
            % read transducer geometry
            h.geom=h5read(filename,[location '/transducer/geom']);
            
            % write time vector <-- to be removed!
            dt=1/h.sampling_frequency;
            h.time=(0:dt:((size(h.data,1)-1)*dt)).'+h.initial_time;
        end
    end
    
    %% implementation methods, to be overloaded
    methods 
        function sig=ir_simple_matlab(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Simple matlab: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_optimized_matlab(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Optimazed matlab: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_mex(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Mex: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_mex_gpu(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Mex gpu: The implementation is not available yet!');
            sig=[];
        end
    end
    
    %% set methods
    methods  
        function set.time(h,input_time)
            assert(size(input_time,1)>size(input_time,2), 'The time vector must be a column vector!')
            h.time=input_time;
            h.sampling_frequency=1/mean(diff(input_time));
            h.initial_time=input_time(1);
        end
        function set.data(h,input_data)
            h.channels=size(input_data,2);
            h.firings=size(input_data,3);
            h.frames=size(input_data,4);
            h.data=input_data;
        end
        function set.geom(h,input_geom)
            assert(ndims(input_geom)==2&&size(input_geom,2)==3, 'Wrong probe geometry definition. It should be a three column matrix [x, y, z]')
            h.geom=input_geom;
        end
    end
    
    %% estimate center frequency
    methods
        function estimate_center_frequency(h)
            switch(h.format)
                case E.signal_format.RF
                    h.center_frequency = tools.estimate_frequency(h.time,h.data);
                case E.signal_format.IQ
                    warning('Cannot estimate center frequency from IQ data, do you want to use the modulation frequency?')
                otherwise
                    error('Unknown signal format!');
            end
        end
    end
end

