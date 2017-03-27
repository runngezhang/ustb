classdef linear_3D_scan < huff.scan
%linear_3D_scan Class defining a linear_3D_scan area. Child of SCAN class 
%
%   See also linear_3D_scan/linear_3D_scan, SCAN/SCAN

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/27$

    properties  (SetAccess = public)
        radial_axis         % Vector containing the coordinates in the radial direction axis [m]
        axial_axis          % Vector containing the coordinates in the axial direction axis [m]
        roll                % Angle between the radial axis and the x-axis [rad]
    end
    
    properties  (Dependent)
        N_radial_axis             % number of pixels in the x_axis
        N_axial_axis              % number of pixels in the z_axis
    end
    
    %% Constructor
    methods (Access = public)
        function h = linear_3D_scan(in_radial_axis,in_axial_axis,in_roll)
            %linear_3D_scan   Constructor of linear_3D_scan class
            %
            %   Syntax:
            %   h = linear_3D_scan(in_radial_axis,in_axial_axis,in_roll)
            %       radial_axis         Vector containing the coordinates in the radial direction axis [m]
            %       axial_axis          Vector containing the coordinates in the axial direction axis [m]
            %       roll                Angle between the radial axis and the x-axis [rad]
            %
            %   See also SCAN
            if nargin>0
                h.radial_axis=in_radial_axis;
            end
            if nargin>1
                h.axial_axis=in_axial_axis;
            end
            if nargin>2
                h.roll=in_roll;
            else
                h.roll=0;
            end
            
            h=h.update_pixel_position();
        end
    end
    
    %% update pixel position
    methods 
        function h=update_pixel_position(h)
            % defining the pixel mesh 
            [R A]=meshgrid(h.radial_axis,h.axial_axis);
            
            % position of the pixels
            if ~isempty(R)&&~isempty(h.roll)
                h.x=R(:)*cos(h.roll);
                h.y=R(:)*sin(h.roll);
                h.z=A(:);
            end
        end
    end
    
    %% Set methods
    methods
        function h=set.radial_axis(h,in_radial_axis)
            assert(size(in_radial_axis,2)==1, 'The input must be a column vector.')
            h.radial_axis=in_radial_axis;
            h=h.update_pixel_position();
        end
        function h=set.axial_axis(h,in_axial_axis)
            assert(size(in_axial_axis,2)==1, 'The input vector must be a column vector.')
            h.axial_axis=in_axial_axis;
            h=h.update_pixel_position();
        end
        function h=set.roll(h,in_roll)
            assert(numel(in_roll)==1, 'The input vector must be a scalar.')
            h.roll=in_roll;
            h=h.update_pixel_position();
        end
    end
    
    %% Get methods
    methods
        function value=get.N_radial_axis(h)
            value=numel(h.radial_axis);
        end
        function value=get.N_axial_axis(h)
            value=numel(h.axial_axis);
        end        
    end
    
end

