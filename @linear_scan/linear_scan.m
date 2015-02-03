classdef linear_scan 
%LINEAR_SCAN Class defining a linear scan area  
%
%   See also RECONSTRUCTION

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $

    properties  (SetAccess = public)
        x_axis      % Vector defining the x coordinates of each row of pixels
        z_axis      % Vector defining the z coordinates of each column of pixels
    end
    
    properties  (SetAccess = protected)
        x_matrix    % Matrix containing the x coordinate of each pixel in the matrix
        z_matrix    % Matrix containing the z coordinate of each pixel in the matrix
        x           % Vector containing the x coordinate of each pixel in the matrix
        z           % Vector containing the z coordinate of each pixel in the matrix
        pixels      % total number of pixels in the matrix
    end
    
    %% Constructor
    methods (Access = public)
        function h = linear_scan(input_x,input_z)
            %LINEAR_SCAN   Constructor of linear_scan class
            %
            %Usage:
            %   h = linear_class(x_vector,z_vector)
            %       x_vector    Vector defining the x coordinates of each row of pixels
            %       z_vector    Vector defining the z coordinates of each column of pixels
            if nargin>0
                h.x_axis=input_x;
            end
            if nargin>1
                h.z_axis=input_z;
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
            
            % scan type
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'linear_scan', 0); 
                H5T.enum_insert (filetype, 'sector_scan', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'scan_type', filetype, space, 'H5P_DEFAULT');
            H5A.write (attr, filetype, uint32(0));          % <---  LINEAR SCAN
                    
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
                        
            % add x-axis
            dset_details.Location = location;
            dset_details.Name = 'x_axis';
            hdf5write(filename, dset_details, h.x_axis, 'WriteMode', 'append');

            % add z-axis
            dset_details.Location = location;
            dset_details.Name = 'z_axis';
            hdf5write(filename, dset_details, h.z_axis, 'WriteMode', 'append');
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
            assert(strcmp(scan_type,'linear_scan'),'Scan format does not match!');

            % x_axis
            h.x_axis=h5read(filename,[location '/x_axis']);
            
            % z_axis
            h.z_axis=h5read(filename,[location '/z_axis']);            
        end
    end
    
    %% Set methods
    methods
        function h=set.x_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The x vector must be a column vector!')
            h.x_axis=input_vector;
            [h.x_matrix, h.z_matrix]=meshgrid(h.x_axis,h.z_axis); 
            h.x=h.x_matrix(:);
            h.z=h.z_matrix(:);
            h.pixels=length(h.x);
        end
        function h=set.z_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The z vector must be a column vector!')
            h.z_axis=input_vector;
            [h.x_matrix, h.z_matrix]=meshgrid(h.x_axis,h.z_axis); 
            h.x=h.x_matrix(:);
            h.z=h.z_matrix(:);
            h.pixels=length(h.x);
        end
    end
end

