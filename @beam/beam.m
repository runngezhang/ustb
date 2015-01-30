classdef beam < handle
%beam   Beam definition
%
%   See also BEAM.BEAM

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $

    properties  (SetAccess = public)
        f_number=1                              % F-number
        apodization=E.apodization_type.none     % Apodization type as given by E.apodization_type
        steer_angle=0                           % Steering angle of the beam (rad)
        smoothing=0                             % Length of the edge smoothing function (0 = no edge smoothing)
        smoothing_order=0                       % Order of the edge smoothing function 
    end
    
    methods (Access = public)
        function h = beam(input_f_number,input_apo,input_steer_angle,input_smooth,input_order)
            %BEAM   Constructor of beam class
            %
            %Usage:
            %   h = beam(f_number,apodization,steer_angle,smoothing,amoothing_order)
            %       f_number=1                              F-number
            %       apodization=E.apodization_type.none     Apodization type as given by E.apodization_type
            %       steer_angle=0                           Steering angle of the beam (rad)
            %       smoothing=0                             Length of the edge smoothing function (0 = no edge smoothing)
            %       smoothing_order=0                       Order of the edge smoothing function 
            
            if nargin>0
                h.f_number=input_f_number;
            end
            if nargin>1
                h.apodization=input_apo;
            end
            if nargin>2
                h.steer_angle=input_steer_angle;
            end
            if nargin>3
                h.smoothing=input_smooth;
            end
            if nargin>4
                h.smoothing_order=input_order;
            end
        end
    end
end

