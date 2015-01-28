classdef vs < dataset
    % Class containig properties and methods to work with Virtual
    % Source (vs) datasets. 
    % vs(name,format,c0,time,data,geom,modulation_frequency) 
            %  * name:      name of the dataset
            %  * format:    format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %  * c0:        reference speed of sound (m/s)
            %  * source:    Position of the virtual sources [x, y, z] (m)
            %  * time:      fast time vector (s)
            %  * data:      Numerical data [time_samples, channels, firings, frames]
            %  * geom:      Probe geometry [x, y, z] (m)
            %  * modulation_frequency:      Modulation frequency (Hz), only for IQ format
    
    properties (SetAccess = public)
        source       % vector containing the angles (rad)
    end
    
    properties (SetAccess = private)
        N           % number of virtual sources in the sequence        
    end
    
    methods (Access = public)        
        function h = vs(name,input_format,input_c0,input_source,input_time,input_data,input_geom,input_modulation_frequency)
            % vs(name,format,c0,source,time,data,geom,modulation_frequency) 
            %  * name:                  Name of the dataset
            %  * format:                Format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %  * c0:                    Reference speed of sound (m/s)
            %  * source:                Position of the virtual sources [x, y, z] (m)
            %  * time:                  Fast time vector (s)
            %  * data:                  Numerical data [time_samples, channels, firings, frames]
            %  * geom:                  Probe geometry [x, y, z] (m)
            %  * modulation_frequency:  Modulation frequency (Hz), only for IQ format
            
            % Call superclass constructor
            h@dataset(name); 
            
            % required data
            h.format = input_format;
            h.c0 = input_c0;
            h.time = input_time;
            h.data = input_data;
            h.geom = input_geom;
            h.source = input_source;
            
            if(h.format==E.signal_format.IQ) 
                h.modulation_frequency=input_modulation_frequency;
            end
            
            % dependent data
            h.Fs = 1/mean(diff(h.time));    % sampling frequency (Hz)
            h.F=size(h.data,4);             % short for the number of frames
            h.N = size(h.source,1);          % short for the number of firing
            h.M = size(h.geom,1);           % short for the number of channels
            h.t0= h.time(1);                % initial time (s)
        end
        
        function [apo]= virtual_apodization(h,beam,r)
            % virtual apodization computed via Alfonso's equation for 
            % mimicking STA transmit apodization
            apo=zeros(h.Nz,h.Nx,h.N);
            Aperture=r.z./beam.FN;
            for na=1:h.N
                xT=r.x-r.z.*(h.source(na,1)-r.x)./(h.source(na,3)-r.z);                            
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
            h.tx_apodization = h.virtual_apodization(beam.transmit,r);
            
            % launch selected implementation
            sig=h.launch_implementation(r,imp);
            
            h.elapsed_time=toc; % elapsed time -> performance test
            
            % image processing
            im=h.envelope(sig);
            %h.show(r,im,60);
        end
    end
    
end

