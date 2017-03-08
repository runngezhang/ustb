classdef linear_scan < scan
%LINEAR_SCAN Class defining a linear_scan area. Child of SCAN class 
%
%   See also LINEAR_SCAN/LINEAR_SCAN, SCAN/SCAN

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/28 $

    properties  (SetAccess = public)
        x_axis           % Vector containing the x coordinates of the x - axis [m]
        z_axis           % Vector containing the z coordinates of the z - axis [m]
    end
    
    properties  (Dependent)
        N_x_axis              % number of pixels in the x_axis
        N_z_axis              % number of pixels in the z_axis
    end
    
    %% Constructor
    methods (Access = public)
        function h = linear_scan(in_x_axis,in_z_axis)
            %LINEAR_SCAN   Constructor of LINEAR_SCAN class
            %
            %   Syntax:
            %   h = linear_scan(x_axis,y_axis)
            %       x    Vector with the x coordinates of each pixel
            %       z    Vector with the z coordinates of each pixel
            %
            %   See also SCAN
            if nargin>0
                h.x_axis=in_x_axis;
            end
            if nargin>1
                h.z_axis=in_z_axis;
            end
            
            h=h.update_pixel_position();
        end
    end
    
    %% update pixel position
    methods 
        function h=update_pixel_position(h)
            % defining the pixel mesh 
            [X Z]=meshgrid(h.x_axis,h.z_axis);
            
            % position of the pixels
            h.x=X(:);
            h.y=0.*X(:);
            h.z=Z(:);
        end
    end
    
    %% Set methods
    methods
        function h=set.x_axis(h,in_x_axis)
            assert(size(in_x_axis,2)==1, 'The input must be a column vector.')
            h.x_axis=in_x_axis;
            h=h.update_pixel_position();
        end
        function h=set.z_axis(h,in_z_axis)
            assert(size(in_z_axis,2)==1, 'The input vector must be a column vector.')
            h.z_axis=in_z_axis;
            h=h.update_pixel_position();
        end
    end
    
    %% Get methods
    methods
        function value=get.N_x_axis(h)
            value=numel(h.x_axis);
        end
        function value=get.N_z_axis(h)
            value=numel(h.z_axis);
        end        
    end
    
%     %% inherit plot method
%     methods
%         function figure_handle=plot(h,figure_handle_in,title_in)
%             figure_handle=plot@scan(h);
%         end
%     end
    
end

