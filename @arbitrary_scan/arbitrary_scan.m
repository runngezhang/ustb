classdef arbitrary_scan 
%ARBITRARY_SCAN Class defining a linear scan area  
%
%   See also RECONSTRUCTION

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/12/07 $

    properties  (SetAccess = public)
        x_matrix           % Matrix containing the x coordinates of each pixel in the matrix
        z_matrix           % Matrix containing the z coordinates of each pixel in the matrix
    end
    
    properties  (SetAccess = private)
        x                  % Vector containing the x coordinate of each pixel in the matrix
        z                  % Vector containing the z coordinate of each pixel in the matrix
        pixels             % total number of pixels in the matrix    
    end
    
    %% Constructor
    methods (Access = public)
        function h = arbitrary_scan(input_x_matrix,input_z_matrix)
            %ARBITRARY_SCAN   Constructor of arbitrary_scan class
            %
            %   Syntax:
            %   h = arbitrary_scan(x,z)
            %       x    Vector with the x coordinates of each pixel
            %       z    Vector with the z coordinates of each pixel
            %
            %   See also ARBITRARY_SCAN
            if nargin>0
                h.x_matrix=input_x_matrix;
            end
            if nargin>1
                h.z_matrix=input_z_matrix;
            end
        end
    end
    
    %% lateral_distance
    methods (Access = public)
        function xd = lateral_distance(h,x0,z0,steer_angle)
            %LATERAL_DISTANCE   Calculates the lateral distance from the center of
            %   the apodization window for a specific scanning mode
            %
            %   Syntax:
            %   h = lateral_distance(element_position,steering_angle)
            %       x0                  Vector containing the x coordinates of the probe elements (either real or virtual) [m]
            %       z0                  Vector containing the x coordinates of the probe elements (either real or virtual) [m]
            %       steering_angle      Steerin angle [rad]
            %
            %   See also LINEAR_SCAN
            xd=abs(x0-h.x+h.z*tan(steer_angle));
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
            %       location                     Name of the group destination
            %
            %   See also RECONSTRUCTION
            
            % scan type
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'linear_scan', 0); 
                H5T.enum_insert (filetype, 'sector_scan', 1); 
                H5T.enum_insert (filetype, 'arbitrary_scan', 2); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'scan_type', filetype, space, 'H5P_DEFAULT');
            H5A.write (attr, filetype, uint32(2));          % <---  ARBITRARY SCAN
                    
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
                        
            % add x-axis
            dset_details.Location = location;
            dset_details.Name = 'x';
            hdf5write(filename, dset_details, h.x, 'WriteMode', 'append');

            % add z-axis
            dset_details.Location = location;
            dset_details.Name = 'z';
            hdf5write(filename, dset_details, h.z, 'WriteMode', 'append');
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
            scan_type=h5readatt(filename,location,'scan_type');
            assert(strcmp(scan_type,'arbitrary_scan'),'Scan format does not match!');

            % x_axis
            h.x=h5read(filename,[location '/x']);
            
            % z_axis
            h.z=h5read(filename,[location '/z']);            
        end
    end
    
    %% Set methods
    methods
        function h=set.x_matrix(h,input_matrix)
            %assert(size(input_vector,1)>size(input_vector,2), 'The x vector must be a column vector!')
            h.x_matrix=input_matrix;
            h.x=h.x_matrix(:);
            h.pixels=min([length(h.x) length(h.z)]);
        end
        function h=set.z_matrix(h,input_matrix)
            %assert(size(input_matrix,1)>size(input_matrix,2), 'The z vector must be a column vector!')
            h.z_matrix=input_matrix;
            h.z=h.z_matrix(:);
            h.pixels=min([length(h.x) length(h.z)]);
        end
    end
end

