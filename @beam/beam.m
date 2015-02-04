classdef beam 
%beam   Beam definition
%
%   See also BEAM.BEAM

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $

    properties  (SetAccess = public)
        f_number=1                              % F-number
        apodization=E.apodization_type.none     % Apodization type as given by E.apodization_type
        steer_angle=0                           % Steering angle of the beam (rad)
        smoothing=0                             % Length of the edge smoothing function (0 = no edge smoothing)
    end
    
    %% constructor
    methods (Access = public)
        function h=beam(input_f_number,input_apo,input_steer_angle,input_smooth)
            %BEAM   Constructor of beam class
            %
            %   Syntax:
            %   h = beam(f_number,apodization,steer_angle,smoothing,amoothing_order)
            %       f_number=1                              F-number
            %       apodization=E.apodization_type.none     Apodization type as given by E.apodization_type
            %       steer_angle=0                           Steering angle of the beam (rad)
            %       smoothing=0                             Length of the edge smoothing function (0 = no edge smoothing)
            %
            %   See also BEAM
            
            if nargin>0
                h.f_number=input_f_number;
            end
            if nargin>1
                h.apodization=input_apo;
            end
            if nargin>2
                h.steer_angle=input_steer_angle;
            end
            if nargin>3
                h.smoothing=input_smooth;
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
            %   See also RECONSTRUCTION
            
            % apodization type
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'none', 0); 
                H5T.enum_insert (filetype, 'boxcar', 1); 
                H5T.enum_insert (filetype, 'hanning', 2); 
                H5T.enum_insert (filetype, 'tukey25', 3); 
                H5T.enum_insert (filetype, 'tukey50', 4); 
                H5T.enum_insert (filetype, 'tukey80', 5); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'apodization', filetype, space, 'H5P_DEFAULT');
            switch (h.apodization)
                case E.apodization_type.none
                    H5A.write (attr, filetype, uint32(0));  
                case E.apodization_type.boxcar
                    H5A.write (attr, filetype, uint32(1));  
                case E.apodization_type.hanning
                    H5A.write (attr, filetype, uint32(2));  
                case E.apodization_type.tukey25
                    H5A.write (attr, filetype, uint32(3));  
                case E.apodization_type.tukey50
                    H5A.write (attr, filetype, uint32(4));  
                case E.apodization_type.tukey80
                    H5A.write (attr, filetype, uint32(5));  
                otherwise
                    error('Unknown apodization type');
            end
                    
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
                        
            % write f_number
            attr = h.f_number;
            attr_details.Name = 'f_number';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');
            
            % write steer angle
            attr = h.steer_angle;
            attr_details.Name = 'steer_angle';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');
            
            % write smoothing
            attr = h.smoothing;
            attr_details.Name = 'smoothing';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');
        end
        
        function h=huff_read(h,filename,location)
            %HUFF_READ    Dumps all the information of the ultrasound
            %dataset to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_READ(file_name,location)
            %       file_name                    Name of the hdf5 file
            %       location                   Name of the group destination
            %
            %   See also US_DATASET
            
            % read signal format 
            apodization=h5readatt(filename,location,'apodization');
            switch(apodization{1})
                case 'none'
                    h.apodization=E.apodization_type.none;
                case 'boxcar'
                    h.apodization=E.apodization_type.boxcar;
                case 'hanning'
                    h.apodization=E.apodization_type.hanning;
                case 'tukey25'
                    h.apodization=E.apodization_type.tukey25;                    
                case 'tukey50'
                    h.apodization=E.apodization_type.tukey50;                    
                case 'tukey80'
                    h.apodization=E.apodization_type.tukey80;                    
                otherwise
                    error('Unknown apodization type!');
            end

            % read f number
            h.f_number=h5readatt(filename,location,'f_number');

            % read angle
            h.steer_angle=h5readatt(filename,location,'steer_angle');

            % read smoothing
            h.smoothing=h5readatt(filename,location,'smoothing');
        end
    end
end

