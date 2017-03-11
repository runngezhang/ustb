classdef matrix_array < probe 
%MATRIX_ARRAY   2D matrix array definition. Child of PROBE class
%
%   See also PROBE, BEAM, PHANTOM

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2017/03/11 $

    %% public properties
    properties  (SetAccess = public)
        pitch_x        % distance between the elements in the azimuth direction [m]
        pitch_y        % distance between the elements in the elevation direction [m]
        N_x            % number of elements in the azimuth direction
        N_y            % number of elements in the elevation direction
        element_width  % width of the elements in the azimuth direction [m]
        element_height % height of the elements in the elevation direction [m]
    end

    %% update method
    methods 
        function h=update(h)
            if ~isempty(h.pitch_x)&&~isempty(h.pitch_y)&&~isempty(h.N_x)&&~isempty(h.N_y) 
                
                if isempty(h.element_width)
                    h.element_width=h.pitch;
                end
                
                if isempty(h.element_height)
                    h.element_height=h.element_width;
                end
                
                % compute element center location
                x0=(1:h.N_x)*h.pitch_x;
                y0=(1:h.N_y)*h.pitch_y;
                x0=x0-mean(x0);
                y0=y0-mean(y0);
                
                % meshgrid
                [X Y]=meshgrid(x,y);

                % assign geometry
                h.geometry=[X(:) Y(:) zeros(h.N_x*h.N_y,3) h.element_width*ones(h.N_x*h.N_y,1) h.element_height*ones(h.N_x*h.N_y,1)]; % probe geometry
            end
        end
    end
    
    %% set methods
    methods  
        function h=set.pitch_x(h,in_pitch)
            assert(numel(in_pitch)==1, 'The input should be a scalar in [m]');
            h.pitch_x=in_pitch;
            h.update();
        end
        function h=set.pitch_y(h,in_pitch)
            assert(numel(in_pitch)==1, 'The input should be a scalar in [m]');
            h.pitch_y=in_pitch;
            h.update();
        end
        function h=set.N_x(h,in_N_elements)
            assert(numel(in_N_elements)==1, 'The input should be a scalar');
            h.N_elements_x=in_N_elements;
            h.update();
        end
        function h=set.N_y(h,in_N_elements)
            assert(numel(in_N_elements)==1, 'The input should be a scalar');
            h.N_elements_y=in_N_elements;
            h.update();
        end
        function h=set.element_width(h,in_width)
            assert(numel(in_width)==1, 'The input should be a scalar in [m]');
            h.element_width=in_width;
            h.update();
        end
        function h=set.element_height(h,in_height)
            assert(numel(in_height)==1, 'The input should be a scalar in [m]');
            h.element_height=in_height;
            h.update();
        end
    end
    
end