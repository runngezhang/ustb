classdef cpw < dataset
    % cpw Class containig properties and methods to work with Coherent
    % plane wave (cpw) datasets.
    
    properties (SetAccess = public)
        angle       % vector containing the angles (rad)
    end
    
    properties (SetAccess = private)
        N           % number of plane waves in the sequence
    end
    
    methods  % formating methods
        function set.angle(h,input_angle)
            if(size(input_angle,1)<size(input_angle,2))
                error('The angle vector must be a column vector!');
            else
                h.angle=input_angle;
            end
        end
    end
    
    methods  (Access = public)       
        function h = cpw(name,input_format,input_c0,input_angle,input_time,input_data,input_geom,input_modulation_frequency)
            % The constructor requires:
            %  * name:      name of the dataset
            %  * format:    format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %  * c0:        reference speed of sound (m/s)
            %  * angle:     vector with the angles in the sequence (rad)
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
            h.angle = input_angle;
            
            if(h.format==E.signal_format.IQ) 
                h.modulation_frequency=input_modulation_frequency;
            end
            
            % dependent data
            h.Fs = 1/mean(diff(h.time));    % sampling frequency (Hz)
            h.F=size(h.data,4);             % short for the number of frames
            h.N = size(h.angle,1);          % short for the number of firings
            h.M = size(h.geom,1);           % short for the number of channels
            h.t0= h.time(1);                % initial time (s)
        end
        
       function [apo]= angular_apodization(h,beam,r)
            % angular apodization computed via Alfonso's equation for 
            % mimicking STA transmit apodization
            apo=zeros(h.Nz,h.Nx,h.N);
            Aperture=r.z./beam.FN;
            for na=1:h.N
                xT=r.x-r.z*tan(h.angle(na));                            
                xd=abs(xT-r.x+r.z*tan(beam.steer_angle)); 
                valid_apo=(xT>h.geom(1,1)).*(xT<h.geom(end,1));
                switch(beam.apodization)
                    case E.apodization_type.none
                        apo(:,:,na)=valid_apo;
                    case E.apodization_type.boxcar
                        apo(:,:,na)=valid_apo.*double(xd<Aperture/2); 
                    case E.apodization_type.hanning
                        apo(:,:,na)=valid_apo.*double(xd<Aperture/2).*(0.5 + 0.5*cos(2*pi*xd./Aperture)); 
                    case E.apodization_type.tukey25
                        roll=0.25;
                        apo(:,:,na)=valid_apo.*(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd/Aperture-roll/2-1/2)));                               
                    case E.apodization_type.tukey50
                        roll=0.5;
                        apo(:,:,na)=valid_apo.*(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd/Aperture-roll/2-1/2)));                               
                    case E.apodization_type.tukey80
                        roll=0.8;
                        apo(:,:,na)=valid_apo.*(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd/Aperture-roll/2-1/2)));                               
                    otherwise
                        error('Unknown apodization type!');
                end
            end
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
            h.tx_apodization = h.angular_apodization(beam.transmit,r);
            
            % launch selected implementation
            sig=h.launch_implementation(r,imp);
            
            h.elapsed_time=toc; % elapsed time -> performance test
            
            % image processing
            im=h.envelope(sig);
            %h.show(r,im,60);
        end
    end
    
end
