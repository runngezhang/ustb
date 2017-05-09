classdef probe 
%probe   Probe definition
%
%   See also BEAM, PHANTOM

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2016/09/01 $

    %% public properties
    properties  (SetAccess = public)
        geometry         % matrix of point scaterers [x y z theta phi width height] - [m m m rad rad m m]
    end
    
    %% dependent properties
    properties  (Dependent)   
        N_elements         % number of elements 
        x                  % center of the element in the x axis[m]
        y                  % center of the element in the y axis[m]
        z                  % center of the element in the z axis[m]
        theta              % orientation of the element in the azimuth direction [rad]
        phi                % orientation of the element in the elevation direction [rad]
        width              % element width [m]
        height             % element height [m]
        r                  % distance from the element center to the origin of coordinates [m]
    end
    
    % UFF variables
    properties (SetAccess = private)
        uff_version = 'v1.0.0';
    end
    
    %% constructor
    methods (Access = public)
        function h=probe(in_geometry)
            %PROBE   Constructor of PROBE class
            %
            %   Syntax:
            %   h = probe(geometry)
            %        geometry         % matrix of point scaterers [x y z theta phi width height] - [m m m rad rad m m]
            %
            %   See also BEAM, PHANTOM
                        
            if nargin>0
               h.geometry=in_geometry;
            end
        end
    end
    
    %% plot methods
    methods
        function figure_handle=plot(h,figure_handle_in,title_in)
            x = [(h.x-h.width/2.*cos(h.theta)).'; (h.x+h.width/2.*cos(h.theta)).'; (h.x+h.width/2.*cos(h.theta)).'; (h.x-h.width/2.*cos(h.theta)).'];
            y = [(h.y-h.height/2.*cos(h.phi)).'; (h.y-h.height/2.*cos(h.phi)).'; (h.y+h.height/2.*cos(h.phi)).'; (h.y+h.height/2.*cos(h.phi)).'; ];
            z = [(h.z+h.width/2.*sin(h.theta)+h.height/2.*sin(h.phi)).'; (h.z-h.width/2.*sin(h.theta)+h.height/2.*sin(h.phi)).'; (h.z-h.width/2.*sin(h.theta)-h.height/2.*sin(h.phi)).'; (h.z+h.width/2.*sin(h.theta)-h.height/2.*sin(h.phi)).'];
            c = linspace(0,1,h.N_elements);
            
            % plotting probe
            if (nargin>1) && ~isempty(figure_handle_in)
                figure_handle=figure(figure_handle_in); hold on;
            else
                figure_handle=figure();
                title('Probe');
            end

            fill3(x*1e3,y*1e3,z*1e3,c); grid on; axis equal tight; hold on;
            plot3(h.x*1e3,h.y*1e3,h.z*1e3,'k+');
            xlabel('x[mm]'); ylabel('y[mm]'); zlabel('z[mm]');
            set(gca,'ZDir','Reverse');
            set(gca,'fontsize',14);
            
            if nargin>2
                title(title_in);
            end
        end
    end
    
    %% set methods
    methods  
        function h=set.geometry(h,in_geometry)
            
            % if element tilt & dimensions aren't provided we set them to 0
            N_columns=size(in_geometry,2);
            if(N_columns<7)
                in_geometry=[in_geometry zeros(size(in_geometry,1),7-size(in_geometry,2))];
            end
            
            assert(size(in_geometry,2)==7, 'The elements matrix should be [x y z theta phi width height] - [m m m rad rad m m]');
            h.geometry=in_geometry;
        end
    end
    
    %% get methods
    methods  
        function value=get.N_elements(h)
            value=size(h.geometry,1);
        end
        function value=get.x(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=h.geometry(:,1);
        end
        function value=get.y(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=h.geometry(:,2);
        end
        function value=get.z(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=h.geometry(:,3);
        end
        function value=get.theta(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=h.geometry(:,4);
        end
        function value=get.phi(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=h.geometry(:,5);
        end
        function value=get.width(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=h.geometry(:,6);
        end
        function value=get.height(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=h.geometry(:,7);
        end
        function value=get.r(h)
            if h.N_elements==0 error('The PROBE class is empty'); end
            value=sqrt(sum(h.geometry(:,1:3).^2,2));
        end
    end
    
     %% Ultrasound File Format
    methods (Access = public)    
        function write(h,filename,location)
            
            %-- Writes all the information to Ultrasound File Format (UFF)
            %-- Syntax:
            %-- write_file_hdf5(file_name)
            %-- file_name: Name of the hdf5 file
            
            %-- write UFF version in the root group
            attr_details.Name = 'version';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.uff_version);
            
            %-- We create a new channel data
            group_name = 'probe';
            fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            gid = H5G.open(fid,[location '/' group_name]);
            
            s_gid = H5G.create(gid,group_name,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
            H5G.close(s_gid);
            H5G.close(gid);
            H5F.close(fid);
            location = location '/' group_name;
   
            uff.append(filename, location, 'geometry', h.sound_speed);
            uff.append(filename, location, 'initial_time', h.initial_time);
            uff.append(filename, location, 'sampling_frequency', h.sampling_frequency);
            uff.append(filename, location, 'modulation_frequency', h.modulation_frequency);
            uff.append(filename, location, 'data', h.data);
            uff.append(filename, location, 'probe', h.probe);

        end
    end
end