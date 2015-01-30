classdef linear_scan < handle
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
    
    %% Protected methods
    methods (Access = protected) 
        function update(h)
            %UPDATE   Recalculates dependent data. Protected method.
            [h.x_matrix, h.z_matrix]=meshgrid(h.x_axis,h.z_axis); 
            h.x=h.x_matrix(:);
            h.z=h.z_matrix(:);
            h.pixels=length(h.x);
        end
    end
    
    %% Set methods
    methods
        function set.x_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The x vector must be a column vector!')
            h.x_axis=input_vector;
            h.update();
        end
        function set.z_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The z vector must be a column vector!')
            h.z_axis=input_vector;
            h.update();
        end
    end
end

