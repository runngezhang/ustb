classdef sector_scan 
%LINEAR_SCAN Class defining a sector scan area  
%
%   See also RECONSTRUCTION

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/02/03 $

    properties  (SetAccess = public)
        apex            % location of the apex from which the sector is defined
        azimuth_axis    % Vector defining the azimuth coordinates of each row of pixels
        depth_axis      % Vector defining the depth coordinates of each column of pixels
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
        function h = sector_scan(azimuth_input,depth_input)
            %SECTOR_SCAN   Constructor of the sector_scan class
            %
            %   Syntax:
            %   h = sector_scan(azimuth_vector,depth_vector)
            %       azimuth_vector    Vector defining the azimuth coordinates of each row of pixels
            %       depth_vector      Vector defining the depth coordinates of each column of pixels
            %
            %   See also SECTOR_SCAN
            
            if nargin>0
                h.azimuth_axis=azimuth_input;
            end
            if nargin>1
                h.depth_axis=depth_input;
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
            %       element_position    Vector containing the x coordinates of the probe elements (either real or virtual) [m]
            %       steering_angle      Steering angle [rad]
            %
            %   See also SECTOR_SCAN
            
            % distance between a point (x0,z0) and a line (x1,z1) -- (x2,z2) 
            x1=h.apex(1);
            z1=h.apex(3);
            x2=h.x;
            z2=h.z;
            
            xd=abs((z2-z1).*x0-(x2-x1).*z0+x2.*z1-z2.*x1)./sqrt((z2-z1).^2+(x2-x1).^2);
            
            %natural_angle=atan2(x0-h.apex(1),z0-h.apex(3));            
            %xd=abs(x0-h.x+h.z.*tan(steer_angle+natural_angle));
        end
    end
    
    %% Spatial step in the beam direction
    methods (Access = public)
        function dz = depth_step(h)
            %DEPTH_STEP   Calculates the spatial step for a given scan area
            %
            %   Syntax:
            %   dz = depth_step()
            %       dz    Spatial step in the beam direction
            %
            %   See also LINEAR_SCAN
            dz=mean(diff(h.depth_axis)); % spatial step in the beam direction
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
            H5A.write (attr, filetype, uint32(1));          % <---  SECTOR SCAN
                    
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);

            % add center
            dset_details.Location = location;
            dset_details.Name = 'center';
            hdf5write(filename, dset_details, h.center, 'WriteMode', 'append');
            
            % add azimuth axis
            dset_details.Location = location;
            dset_details.Name = 'azimuth_axis';
            hdf5write(filename, dset_details, h.x_axis, 'WriteMode', 'append');

            % add depth axis
            dset_details.Location = location;
            dset_details.Name = 'depth_axis';
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
            assert(strcmp(scan_type,'sector_scan'),'Scan format does not match!');

            % azimuth axis
            h.center=h5read(filename,[location '/center']);
            
            % azimuth axis
            h.x_axis=h5read(filename,[location '/azimuth_axis']);
            
            % depth axis
            h.z_axis=h5read(filename,[location '/depth_axis']);            
        end
    end
    
    %% Set methods
    methods
        function h=set.azimuth_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The azimuth vector must be a column vector!')
            h.azimuth_axis=input_vector;
            
            if ~isempty(h.azimuth_axis)
                [aa, dd]=meshgrid(h.azimuth_axis,h.depth_axis); 
            
                h.x_matrix=dd.*sin(aa)+h.apex(1);
                h.z_matrix=dd.*cos(aa)+h.apex(3);

                h.x=h.x_matrix(:);
                h.z=h.z_matrix(:);
                h.pixels=length(h.x);
            end
        end
        function h=set.depth_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The depth vector must be a column vector!')
            h.depth_axis=input_vector;
            
            if ~isempty(h.azimuth_axis)
                [aa, dd]=meshgrid(h.azimuth_axis,h.depth_axis); 

                h.x_matrix=dd.*sin(aa)+h.apex(1);
                h.z_matrix=dd.*cos(aa)+h.apex(3);

                h.x=h.x_matrix(:);
                h.z=h.z_matrix(:);
                h.pixels=length(h.x);
            end
        end
    end
end

