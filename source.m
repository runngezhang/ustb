classdef source 
%source   Source definition
%

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/21 $

    %% public properties
    properties  (SetAccess = public)
        distance     % distance from the source location to the origin of coordinates [m]
        azimuth      % angle from the source location to the plane YZ [rad]
        elevation    % angle from the source location to the plane XZ [rad]
    end
    
    %% dependent properties
    properties  (Dependent)
        xyz          % location of the source [m m m] if the source is not at infinity        
        x
        y
        z
    end
    
    %% constructor
    methods (Access = public)
        function h=source(distance_in,azimuth_in,elevation_in)
            %SOURCE   Constructor of SOURCE class
            %
            %   Syntax:
            %   h = source(distance,azimuth,elevation)
            %           distance     % distance from the source location to the origin of coordinates [m]
            %           azimuth      % angle from the source location to the plane YZ [rad]
            %           elevation    % angle from the source location to the plane XZ [rad]
            %
            %   All input parameters can be inserted after declaration.
            %
            %   See also SOURCE
            
            if nargin>0
               h.distance=distance_in;
            end
            if nargin>1
               h.azimuth=azimuth_in;
            end       
            if nargin>2
               h.elevation=elevation_in;
            end              
            if nargin>3
               error('Too many input parameters.');
            end              
        end
    end
    
    %% plot methods
    methods
        function plot(h)
%           figure; 
%           plot3(1:numel(h.delay),h.delay*1e6); grid on; axis tight;
%           xlabel('element');
%           ylabel('delay [\mus]');
%           set(gca,'fontsize',14);
%           ylim([min([h.delay*1e6; -1]) max([1; h.delay*1e6])]);
%           title('Delays');
        end
    end
    
    %% set methods
    methods  
         function h=set.distance(h,in_distance)
             assert(numel(in_distance)==1, 'The distance should be a escalar in [m]');
             h.distance=in_distance;
         end
         function h=set.azimuth(h,in_azimuth)
             assert(numel(in_azimuth)==1, 'The azimuth should be a escalar in [rad]');
             h.azimuth=in_azimuth;
         end
         function h=set.elevation(h,in_elevation)
             assert(numel(in_elevation)==1, 'The elevation should be a escalar in [rad]');
             h.elevation=in_elevation;
         end
         function h=set.xyz(h,in_xyz)
             assert(size(in_xyz,2)==3, 'The xyz must be an array [x y z] - [m m m]');
             h.x=in_xyz(1);
             h.y=in_xyz(2);
             h.z=in_xyz(3);
             
             h.distance=norm(h.xyz,2);
             h.azimuth=atan2(h.x,h.z);
             h.elevation=atan2(h.y,h.z);
         end
         function h=set.x(h,in_x)
             assert(numel(in_x)==1, 'The x must be an scalar in [m]');
             h.x=in_x;
             h.distance=norm(h.xyz,2);
             h.azimuth=atan2(h.x,h.z);
             h.elevation=atan2(h.y,h.z);
         end
         function h=set.y(h,in_y)
             assert(numel(in_y)==1, 'The y must be an scalar in [m]');
             h.x=in_y;
             h.distance=norm(h.xyz,2);
             h.azimuth=atan2(h.x,h.z);
             h.elevation=atan2(h.y,h.z);
         end
         function h=set.z(h,in_z)
             assert(numel(in_z)==1, 'The z must be an scalar in [m]');
             h.z=in_z;
             h.distance=norm(h.xyz,2);
             h.azimuth=atan2(h.x,h.z);
             h.elevation=atan2(h.y,h.z);
         end         
    end
    
    %% get methods
    methods  
        function value=get.x(h)
            value=h.distance*cos(h.azimuth)*sin(h.elevation);
        end
        function value=get.y(h)
            value=h.distance*sin(h.azimuth)*sin(h.elevation);
        end
        function value=get.z(h)
            value=h.distance*cos(h.elevation);
        end
        function value=get.xyz(h)
            value=[h.x h.y h.z];
        end
    end
    
    %% plot methods
    methods
        function plot(h)
            % plotting phantom
            if isInf(h.distance)
            else
                figure;
                plot3(h.x*1e3,h.y*1e3,h.z*1e3,'r.'); grid on; axis equal;
                xlabel('x[mm]'); ylabel('y[mm]'); zlabel('z[mm]');
                set(gca,'ZDir','Reverse');
                set(gca,'fontsize',14);
                title('Source');
            end
        end
    end
end