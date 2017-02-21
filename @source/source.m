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
%         function h=set.delay(h,in_delay)
%             assert(size(in_delay,2)==1, 'The delays should be a column vector [N_elements 1]');
%             h.delay=in_delay;
%         end
    end
    
    %% get methods
    methods  
%         function value=get.N_elements(h)
%             value=numel(h.delay);
%         end
    end
    
end