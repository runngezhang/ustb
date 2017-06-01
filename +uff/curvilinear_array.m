classdef curvilinear_array < uff.probe 
%CURVILINEAR_ARRAY   Curvilinear array definition. Child of PROBE class
%
%   See also PROBE, BEAM, PHANTOM

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2017/06/01 $

    %% public properties
    properties  (SetAccess = public)
        pitch          % distance between the elements in the radial direction [m]
        N              % number of elements
        radius         % radius of the curvilinear array [m]
        element_width  % width of the elements in the azimuth direction [m]
        element_height % height of the elements in the elevation direction [m]
    end
    
     properties  (Dependent)
        maximum_angle  % angle of the outermost elements in the array
    end 

    %% update method
    methods 
        function h=update(h)
            if ~isempty(h.pitch)&~isempty(h.N)&~isempty(h.radius) 
                
                if isempty(h.element_width)
                    h.element_width=h.pitch;
                end
                
                if isempty(h.element_height)
                    h.element_height=10*h.element_width;
                end

                % compute element coordinates
                dtheta=2*asin(h.pitch/2/h.radius); 
                theta=(0:h.N-1)*dtheta; theta=theta-mean(theta);
                x0=h.radius*sin(theta);
                z0=h.radius*cos(theta)-h.radius;

                % assign geometry
                h.geometry=[x0(:) zeros(h.N,1) z0(:) theta(:) zeros(h.N,1) h.element_width*ones(h.N,1) h.element_height*ones(h.N,1)]; % probe geometry
            end
        end
    end
    
    %% set methods
    methods  
        function h=set.pitch(h,in_pitch)
            assert(numel(in_pitch)==1, 'The input should be a scalar in [m]');
            h.pitch=in_pitch;
            h=h.update();
        end
        function h=set.N(h,in_N_elements)
            assert(numel(in_N_elements)==1, 'The input should be a scalar');
            h.N=in_N_elements;
            h=h.update();
        end
        function h=set.element_width(h,in_width)
            assert(numel(in_width)==1, 'The input should be a scalar in [m]');
            h.element_width=in_width;
            h=h.update();
        end
        function h=set.element_height(h,in_height)
            assert(numel(in_height)==1, 'The input should be a scalar in [m]');
            h.element_height=in_height;
            h=h.update();
        end
        function h=set.radius(h,in_radius)
            assert(numel(in_radius)==1, 'The input should be a scalar in [m]');
            h.radius=in_radius;
            h=h.update();
        end
    end
     %% Get methods
    methods
        function value=get.maximum_angle(h)
            value=max(abs(h.geometry(:,4)));
        end
    end
end