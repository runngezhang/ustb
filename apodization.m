classdef apodization 
%APDOIZATION   apodization definition
%
%   See also SOURCE, PHANTOM, PROBE, PULSE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/07 $

    %% public properties
    properties  (SetAccess = public)
        probe               % PROBE class
        apex                % SOURCE class
        scan                % SCAN class        
        window              % apodization window
        f_number            % F-number [Fx Fy] [unitless unitless] 
        tilt                % tilt angle [azimuth elevation] [rad rad] 
    end
    
    %% dependent properties
    properties  (Dependent)   
        data                % apodization data 
    end
    
    %% constructor
    methods (Access = public)
        function h=apodization()
            %apodization   Constructor of apodization class
            %
            %   Syntax:
            %   h = apodization()
            %           apex      % apex class
            %
            %   See also apodization, apex, PHANTOM, PROBE, PULSE
            
            h.window=window.none;
            h.probe=probe();
            h.apex=source();
            h.scan=scan(0,0,0);
            h.tilt=[0 0];
        end
    end
    
    %% set methods
    methods  
        function h=set.apex(h,in_source)
            assert(isa(in_source,'source'), 'The input apex is not a SOURCE class. Check HELP apex');
            h.apex=in_source;
        end
        function h=set.probe(h,in_probe)
            assert(isa(in_probe,'probe'), 'The input probe is not a PROBE class. Check HELP PROBE');
            h.probe=in_probe;
        end
        function h=set.scan(h,in_scan)
            assert(isa(in_scan,'scan'), 'The input scan is not a SCAN class. Check HELP SCAN');
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
            assert(size(in_tilt,1)==1,size(in_tilt,1)==2,'The tilt must be a row vector [azimuth elevation]');
            h.tilt=in_tilt;
        end  
        function h=set.window(h,in_window)
             assert(strcmp(class(in_window),'window'),'The window input should be a WINDOW class. Check help WINDOW');
             h.window=in_window;
        end
    end
    
    %% get method
    methods
        function value=get.data(h)
           
            % checking we have all we need
            assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(numel(h.scan)>0,'The SCAN parameter is not set.');
            assert(numel(h.apex)>0,'The APEX parameter is not set.');
            %assert(all(h.scan.z>0),'Cannot compute apodization for points behind the XY plane');
            
            % NONE APODIZATION
            if(h.window==window.none)
                value=ones(h.scan.N_pixels,h.probe.N_elements);
                return;
            end
            
            % STA APODIZATION (just the element closest to apex)
            if (h.window==window.sta)
                dist=sqrt((h.probe.x-h.apex.x).^2+(h.probe.y-h.apex.y).^2+(h.probe.z-h.apex.z).^2);
                value=ones(h.scan.N_pixels,1)*double(dist==min(dist(:)));
                return;
            end
            
            % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
            if isinf(h.apex.distance)
                origin_x=h.scan.x+h.scan.z*tan(h.apex.azimuth)+h.scan.z*tan(h.tilt(1));
                origin_y=h.scan.y+h.scan.z*tan(h.apex.elevation)+h.scan.z*tan(h.tilt(2));
            else
                origin_x=h.apex.x-h.apex.z*(h.scan.x-h.apex.x)./(h.scan.z-h.apex.z)+h.scan.z*tan(h.tilt(1));
                origin_y=(h.scan.y-h.apex.y)./(h.scan.x-h.apex.x).*(origin_x-h.apex.x)+h.scan.y+h.scan.z*tan(h.tilt(2));
            end
            
            origin_x(isnan(origin_x))=0; % solve division by 0
            origin_y(isnan(origin_y))=0; % solve division by 0
            
            x_dist=abs(ones(h.scan.N_pixels,1)*(h.probe.x.') - origin_x*ones(1,h.probe.N_elements));
            y_dist=abs(ones(h.scan.N_pixels,1)*(h.probe.y.') - origin_y*ones(1,h.probe.N_elements));    
            
            % computing aperture
            Aperture_x=(abs(h.scan.z)./h.f_number(1))*ones(1,h.probe.N_elements);
            Aperture_y=(abs(h.scan.z)./h.f_number(2))*ones(1,h.probe.N_elements);

            % SWITCH
            switch(h.window)
                % BOXCAR/FLAT/RECTANGULAR
                case window.boxcar
                    value=double(x_dist<=Aperture_x/2).*...
                          double(y_dist<=Aperture_y/2);
                % HANNING
                case window.hanning
                    value=double(x_dist<=Aperture_x/2).*(0.5 + 0.5*cos(2*pi*x_dist./Aperture_x)).*... 
                          double(y_dist<=Aperture_y/2).*(0.5 + 0.5*cos(2*pi*y_dist./Aperture_y));
                % HAMMING
                case window.hamming
                    value=double(x_dist<=Aperture_x/2).*(0.53836 + 0.46164*cos(2*pi*x_dist./Aperture_x)).*...
                          double(y_dist<=Aperture_y/2).*(0.53836 + 0.46164*cos(2*pi*y_dist./Aperture_y));
                % TUKEY25        
                case window.tukey25
                    roll=0.25;
                    value=(x_dist<=(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<=(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY50        
                case window.tukey50
                    roll=0.50;
                    value=(x_dist<=(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<=(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY75        
                case window.tukey75
                    roll=0.75;
                    value=(x_dist<=(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<=(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY80        
                case window.tukey80
                    roll=0.80;
                    value=(x_dist<=(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<=(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                otherwise
                    error('Unknown apodization type!');
            end
            
%             % edge smoothing
%             if(size(apo,2)>(2*(beam.smoothing+1)))
%                 % compute edge smoothing mask
%                 mask=0.5-0.5*cos((0:beam.smoothing)/beam.smoothing*pi);  % vector mask
%                 mask=(ones(h.scan.pixels,1)*mask(2:(beam.smoothing+1))); % matrix mask
%                 
%                 % modifying apodization edges
%                 apo(:,1:beam.smoothing)=apo(:,1:beam.smoothing).*mask;
%                 apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events)=apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events).*fliplr(mask);
%             end
        end
    end
end