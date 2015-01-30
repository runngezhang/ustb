classdef sta < dataset
%STA    Synthetic transmit aperture (sta) dataset.
%
%   See also DATASET, STA.STA

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $

    properties (SetAccess = public)
    end
    
    properties (SetAccess = private)
    end
    
    %% constructor
    methods (Access = public)        
        function h = sta(name,input_format,input_c0,input_time,input_data,input_geom,input_modulation_frequency)
            %STA    Constructor of the sta class.
            %
            %   Syntax:
            %   STA(name,format,c0,time,data,geom,modulation_frequency) 
            %       name                    Name of the dataset
            %       format                  Format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %       c0                      Reference speed of sound (m/s)
            %       time                    Time vector (s)
            %       data                    Numerical data [time_samples, channels, firings, frames]
            %       geom                    Probe geometry [x, y, z] (m)
            %       modulation_frequency    Modulation frequency (Hz) - required only for IQ format
            %
            %   See also STA, DATASET
                        
            % call superclass constructor
            h@dataset(name); 
            
            % required data
            h.format = input_format;
            h.c0 = input_c0;
            h.time = input_time;
            h.data = input_data;
            h.geom = input_geom;
            if(h.format==E.signal_format.IQ) 
                h.modulation_frequency=input_modulation_frequency;
            end
            
            % checks
            assert(h.channels==size(h.geom,1),'The length of the geometry vector should match number of channels in the dataset');
            assert(h.channels==h.firings,'The number of channels must match the number of firings in STA');
        end
    end
    
    %% image recontructor
    methods (Access = public)        
        function elapsed_time=image_reconstruction(h,recons,implem)
            %IMAGE_RECONSTRUCTION    Method for image reconstruction of
            %Synthetic transmit aperture datasets.
            %
            %   Syntax:
            %   IMAGE_RECONSTRUCTION(reconstruction,implementation)
            %       reconstruction          Class containing the reconstruction specification
            %       implementation          Enumeration specifying the algorithm 
            %
            %   See also RECONSTRUCTION, STA
                        
            % default implementation
            if ~exist('implem') implem=E.implementation.mex; end
            
            % initial time -> performance test
            tic; 
            
            % precompute transmit and receive apodization
            xT=ones(recons.scan.pixels,1)*(h.geom(:,1).');   % position of transmit and receive element
            h.transmit_apodization = recons.calculate_apodization(recons.transmit_beam,xT);
            h.receive_apodization = recons.calculate_apodization(recons.receive_beam,xT);
            
            % launch selected implementation
            recons.data=h.launch_implementation(recons,implem);
            
            % copy data to reconstruction
            recons.format=h.format;
            [recons.central_frequency recons.bandwidth]=tools.estimate_frequency(h.time,h.data);
                        
            % elapsed time for performance test
            elapsed_time=toc; 
        end
    end
    
end

