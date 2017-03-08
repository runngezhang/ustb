classdef scan 
%scan Class defining a scan area  
%
%   See also WAVE, RAW_DATA, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/24 $

    properties  (SetAccess = public)
        x                  % Vector containing the x coordinate of each pixel in the matrix
        y                  % Vector containing the x coordinate of each pixel in the matrix
        z                  % Vector containing the z coordinate of each pixel in the matrix
    end
    
    properties  (Dependent)
        N_pixels           % total number of pixels in the matrix    
    end
    
    %% Constructor
    methods (Access = public)
        function h = scan(in_x,in_y,in_z)
            %scan   Constructor of scan class
            %
            %   Syntax:
            %   h = scan(x,y,z)
            %       x    Vector with the x coordinates of each pixel
            %       z    Vector with the z coordinates of each pixel
            %
            %   See also scan
            if nargin>0
                h.x=in_x;
            end
            if nargin>1
                h.y=in_y;
            end            
            if nargin>2
                h.z=in_z;
            end
        end
    end
    
    %% plot methods
    methods
        function figure_handle=plot(h,figure_handle_in,title_in)
            % plotting scan
            if (nargin>1) && ~isempty(figure_handle_in)
                figure_handle=figure(figure_handle_in); hold on;
            else
                figure_handle=figure();
                title('Probe');
            end

            plot3(h.x*1e3,h.y*1e3,h.z*1e3,'k.');
            xlabel('x[mm]'); ylabel('y[mm]'); zlabel('z[mm]');
            set(gca,'ZDir','Reverse');
            set(gca,'fontsize',14);
            
            if nargin>2
                title(title_in);
            end
        end
    end
    
    %% Set methods
    methods
        function h=set.x(h,in_x)
            assert(size(in_x,2)==1, 'The x vector must be a column vector.')
            h.x=in_x;
        end
        function h=set.y(h,in_y)
            assert(size(in_y,2)==1, 'The y vector must be a column vector.')
            h.y=in_y;
        end
        function h=set.z(h,in_z)
            assert(size(in_z,2)==1, 'The z vector must be a column vector.')
            h.z=in_z;
        end
    end
    
    %% Get methods
    methods
        function value=get.N_pixels(h)
            value=min([numel(h.x) numel(h.y) numel(h.z)]);
        end
    end
end

