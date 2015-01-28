classdef sta < dataset
    % sta Class containig properties and methods to work with Synthetic
    % transmit aperture (sta) dataset.
    
    properties (SetAccess = public)
    end
    
    properties (SetAccess = private)
    end
    
    methods (Access = public)        
        function h = sta(name,input_format,input_c0,input_time,input_data,input_geom,input_modulation_frequency)
            % sta(name,format,c0,time,data,geom,modulation_frequency) 
            %  * name:      name of the dataset
            %  * format:    format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %  * c0:        reference speed of sound (m/s)
            %  * time:      fast time vector (s)
            %  * data:      Numerical data [time_samples, channels, firings, frames]
            %  * geom:      Probe geometry [x, y, z] (m)
            %  * modulation_frequency:      Modulation frequency (Hz), only for IQ format
            
            % Call superclass constructor
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
            
            % dependent data
            h.Fs = 1/mean(diff(h.time));
            h.M = size(h.geom,1);
            h.F=size(h.data,4);  % short for the number of frames
            h.t0= h.time(1);
        end
        
        function [sig,im] = image_reconstruction(h,beam,r,imp)
            % [sig,im] = image_reconstruction(beam,r,imp)
            if (nargin<4) imp=E.implementation.mex; end
            
            % storing this makes the code cleaner
            h.Nz=size(r.x,1);    % number of pixels in the axial direction
            h.Nx=size(r.x,2);    % number of pixels in the lateral direction
            
            tic; % initial time -> performance test
            
            % set tx and rx apodization
            h.rx_apodization = h.linear_apodization(beam.receive,r);
            h.tx_apodization = h.linear_apodization(beam.transmit,r);
            
            % launch selected implementation
            sig=h.launch_implementation(r,imp);
            
            h.elapsed_time=toc; % elapsed time -> performance test
            
            % image processing
            im=h.envelope(sig);
            %h.show(r,im,60);
        end
    end
    
end

