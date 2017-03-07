classdef apodizator 
%APODIZATOR   apodizator definition
%
%   See also SOURCE, PHANTOM, PROBE, PULSE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/07 $

    %% public properties
    properties  (SetAccess = public)
        probe               % PROBE class
        origin              % SOURCE class
        scan                % SCAN class        
        apodization_window  % apodization window
        f_number            % F-number [Fx Fy] [unitless unitless] 
        tilt                % tilt angle [azimuth elevation] [rad rad] 
    end
    
    %% constructor
    methods (Access = public)
        function h=apodizator()
            %apodizator   Constructor of apodizator class
            %
            %   Syntax:
            %   h = apodizator()
            %           origin      % origin class
            %
            %   See also apodizator, origin, PHANTOM, PROBE, PULSE
            
            h.apodization_window=window.none;
            h.probe=probe();
            h.origin=source();
            h.scan=scan(0,0,0);
        end
    end
    
    %% set methods
    methods  
        function h=set.origin(h,in_source)
            assert(strcmp(class(in_source),'source'), 'The input origin is not a SOURCE class. Check HELP origin');
            h.origin=in_source;
        end
        function h=set.probe(h,in_probe)
            assert(strcmp(class(in_probe),'probe'), 'The input probe is not a PROBE class. Check HELP PROBE');
            h.probe=in_probe;
        end
        function h=set.scan(h,in_scan)
            assert(strcmp(class(in_scan),'scan'), 'The input scan is not a SCAN class. Check HELP SCAN');
            h.scan=in_scan;
        end
        function h=set.f_number(h,in_f_number)
            if(numel(in_f_number)==1) % we allow for escalar input
                 in_f_number=[in_f_number in_f_number];
             end 
            assert(size(in_f_number,1)==1&&size(in_f_number,2)==2,'The f-number must be a row vector [Fx Fy]');
            h.f_number=in_f_number;
        end        
        function h=set.tilt(h,in_tilt)
            if(numel(in_tilt)==1) % we allow for escalar input
                 in_tilt=[in_tilt 0];
             end 
            assert(size(in_tilt)==[1 2],'The tilt must be a row vector [azimuth elevation]');
            h.tilt=in_tilt;
        end  
        function h=set.apodization_window(h,in_window)
             assert(strcmp(class(in_window),'window'),'The apodization_window input should be a WINDOW class. Check help WINDOW');
             h.apodization_window=in_window;
        end
    end
    
    %% go method
    methods
        function value=go(h)
            % NONE APODIZATION
            if(h.apodization_window==window.none)
                value=ones(h.scan.N_pixels,h.probe.N_elements);
                return;
            end
            
            % compute lateral distance
            x_dist=abs(h.probe.x-h.origin.x).';    
            y_dist=abs(h.probe.y-h.origin.y).';    

            % STA APODIZATION (just element closest to origin)
            if (h.apodization_window==window.sta)
                value=double(x_dist==min(x_dist)).*double(y_dist==min(y_dist));
                return;
            end
            
            % computing aperture
            Aperture_x=h.scan.z./h.f_number(1);
            Aperture_y=h.scan.z./h.f_number(2);

            % SWITCH
            switch(h.apodization_window)
                % BOXCAR/FLAT/RECTANGULAR
                case window.boxcar
                    value=double(x_dist<Aperture_x/2).*double(y_dist<Aperture_y/2);
                % HANNING
                case window.hanning
                    value=double(x_dist<Aperture_x/2).*(0.5 + 0.5*cos(2*pi*x_dist./Aperture_x)).*... 
                          double(y_dist<Aperture_y/2).*(0.5 + 0.5*cos(2*pi*y_dist./Aperture_y));
                % HAMMING
                case window.hamming
                    value=double(x_dist<Aperture_x/2).*(0.53836 + 0.46164*cos(2*pi*x_dist./Aperture_x)).*...
                          double(y_dist<Aperture_y/2).*(0.53836 + 0.46164*cos(2*pi*y_dist./Aperture_y));
                % TUKEY25        
                case window.tukey25
                    roll=0.25;
                    value=(x_dist<(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY50        
                case window.tukey50
                    roll=0.50;
                    value=(x_dist<(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY75        
                case window.tukey75
                    roll=0.75;
                    value=(x_dist<(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY80        
                case window.tukey80
                    roll=0.80;
                    value=(x_dist<(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                otherwise
                    error('Unknown apodizator type!');
            end
            
%             % edge smoothing
%             if(size(apo,2)>(2*(beam.smoothing+1)))
%                 % compute edge smoothing mask
%                 mask=0.5-0.5*cos((0:beam.smoothing)/beam.smoothing*pi);  % vector mask
%                 mask=(ones(h.scan.pixels,1)*mask(2:(beam.smoothing+1))); % matrix mask
%                 
%                 % modifying apodizator edges
%                 apo(:,1:beam.smoothing)=apo(:,1:beam.smoothing).*mask;
%                 apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events)=apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events).*fliplr(mask);
%             end
            
        end
    end
end