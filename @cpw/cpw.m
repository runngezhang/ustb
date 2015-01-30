classdef cpw < dataset
%CPW    Coherent plane wave (sta) dataset.
%
%   See also DATASET, CPW.CPW

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $
    
    properties (SetAccess = public)
        angle       % vector containing the angles (rad)
    end
    
    %% Constructor
    methods  (Access = public)       
        function h = cpw(name,input_format,input_c0,input_angle,input_time,input_data,input_geom,input_modulation_frequency)
            %CPW    Constructor of the sta class.
            %
            %   Syntax:
            %   CPW(name,format,c0,angle,time,data,geom,modulation_frequency) 
            %       name                    Name of the dataset
            %       format                  Format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %       c0                      Reference speed of sound (m/s)
            %       angle                   Plane wave angle vector (rad)
            %       time                    Time vector (s)
            %       data                    Numerical data [time_samples, channels, firings, frames]
            %       geom                    Probe geometry [x, y, z] (m)
            %       modulation_frequency    Modulation frequency (Hz) - required only for IQ format
            %
            %   See also CPW, DATASET
        
            % Call superclass constructor
            h@dataset(name); 
            
            % required data
            h.format = input_format;
            h.c0 = input_c0;
            h.time = input_time;
            h.data = input_data;
            h.geom = input_geom;
            h.angle = input_angle;
            
            if(h.format==E.signal_format.IQ) 
                h.modulation_frequency=input_modulation_frequency;
            end
            
            % checks
            assert(h.channels==size(h.geom,1),'The length of the geometry vector should match number of channels in the dataset');
            assert(h.firings==length(h.angle),'The angle definition should match the number of firings in the dataset');
        end
    end

    %% image reconstruction
    methods  (Access = public)       
         function elapsed_time = image_reconstruction(h,recons,implem)
            %IMAGE_RECONSTRUCTION    Method for image reconstruction of Coherent plane wave datasets
            %
            %   Syntax:
            %   IMAGE_RECONSTRUCTION(transmit_beam,receive_beam,scan,implementation)
            %       reconstruction          Class containing the reconstruction specification
            %       implementation          Enumeration specifying the algorithm 
            %
            %   See also RECONSTRUCTION, CPW
             
            % default implementation
            if ~exist('implem') implem=E.implementation.mex; end
            
            % initial time -> performance test
            tic; 
            
            % precompute transmit and receive apodization
            xT=recons.scan.x*ones(1,h.firings)-recons.scan.z*tan(h.angle.'); % position of equivalent receive element -> Alfonso's equation 
            valid_apodization =(xT>h.geom(1,1))&(xT<h.geom(end,1));         % check we don't get out of the aperture
            h.transmit_apodization = valid_apodization.*recons.calculate_apodization(recons.transmit_beam,xT);
            xR=ones(recons.scan.pixels,1)*(h.geom(:,1).');                  % position of receive element
            h.receive_apodization = recons.calculate_apodization(recons.receive_beam,xR);
            
            % launch selected implementation
            recons.data=h.launch_implementation(recons,implem);
            
            % copy data to reconstruction
            recons.format=h.format;
            [recons.central_frequency recons.bandwidth]=tools.estimate_frequency(h.time,h.data);
                        
            % elapsed time for performance test
            elapsed_time=toc; 
        end
    end

    %% set methods
    methods  
        function set.angle(h,input_angle)
            assert(size(input_angle,1)>=size(input_angle,2), 'The angle must be a column vector');
            h.angle=input_angle;
        end
    end
    
end
