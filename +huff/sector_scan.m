classdef sector_scan < huff.scan
%SECTOR_SCAN Class defining a SECTOR_SCAN area. Child of SCAN class 
%
%   See also SECTOR_SCAN/SECTOR_SCAN, SCAN/SCAN

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/14 $

    properties  (SetAccess = public)
        azimuth_axis     % Vector containing the azimuth coordinates of the azimuth axis [rad]
        depth_axis       % Vector containing the distance coordinates of the distance axis [m]
        apex             % POINT class
    end
    
    properties  (Dependent)
        N_azimuth_axis            % number of pixels in the x_axis
        N_depth_axis              % number of pixels in the z_axis
    end
    
    %% Constructor
    methods (Access = public)
        function h = sector_scan(in_azimuth_axis,in_depth_axis,apex)
            %sector_scan   Constructor of sector_scan class
            %
            %   Syntax:
            %   h = sector_scan(x_axis,y_axis)
            %       x    Vector with the x coordinates of each pixel
            %       z    Vector with the z coordinates of each pixel
            %
            %   See also SCAN

            % this goes first, otherwise we cannot update
            if nargin>2
                h.apex=in_apex;
            else
                h.apex=huff.point();
                h.apex.xyz=[0 0 0];
            end

            if nargin>0
                h.azimuth_axis=in_azimuth_axis;
            end
            if nargin>1
                h.depth_axis=in_depth_axis;
            end
            
            h=h.update_pixel_position();
        end
    end
    
    %% update pixel position
    methods 
        function h=update_pixel_position(h)
            % defining the pixel mesh 
            [T R]=meshgrid(h.azimuth_axis,h.depth_axis);
            
            % position of the pixels
            h.x=R(:).*sin(T(:))+h.apex.x;
            h.y=0.*R(:)+h.apex.y;
            h.z=R(:).*cos(T(:))+h.apex.z;
        end
    end
    
    %% Set methods
    methods
        function h=set.azimuth_axis(h,in_azimuth_axis)
            assert(size(in_azimuth_axis,2)==1, 'The input must be a column vector.')
            h.azimuth_axis=in_azimuth_axis;
            h=h.update_pixel_position();
        end
        function h=set.depth_axis(h,in_depth_axis)
            assert(size(in_depth_axis,2)==1, 'The input vector must be a column vector.')
            h.depth_axis=in_depth_axis;
            h=h.update_pixel_position();
        end
        function h=set.apex(h,in_apex)
            assert(isa(in_apex,'huff.point'), 'The input is not a SOURCE class. Check HELP SOURCE');
            h.apex=in_apex;
        end
    end
    
    %% Get methods
    methods
        function value=get.N_azimuth_axis(h)
            value=numel(h.azimuth_axis);
        end
        function value=get.N_depth_axis(h)
            value=numel(h.depth_axis);
        end        
    end
   
end

