classdef apodization 
%APDOIZATION   apodization definition
%
%   See also POINT, PHANTOM, PROBE, PULSE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/03/07 $

    %% public properties
    properties  (SetAccess = public)
        probe               % PROBE class
        apex                % POINT class
        scan                % SCAN class        
        window              % apodization window
        f_number            % F-number [Fx Fy] [unitless unitless] 
        sequence            % collection of WAVE classes (for computation of transmit apodization in synthetic scenarios)
        tilt                % tilt angle [azimuth elevation] [rad rad] 
    end
    
    %% dependent properties
    properties  (Dependent)   
        data                        % apodization data 
        origin                      % beam origin
        propagation_distance        % distance from the origin to the pixel
        element_position_matrix     % location of the elements as seen by the pixels [x y z] [m m m]
        N_elements                  % number of elements (real or synthetic)
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
            
            h.window=huff.window.none;
            h.probe=huff.probe();
            h.apex=huff.point();
            h.scan=huff.scan(0,0,0);
            h.tilt=[0 0];
        end
    end
    
    %% set methods
    methods  
        function h=set.apex(h,in_source)
            assert(isa(in_source,'huff.point'), 'The input apex is not a POINT class. Check HELP POINT');
            h.apex=in_source;
        end
        function h=set.probe(h,in_probe)
            assert(isa(in_probe,'huff.probe'), 'The input probe is not a PROBE class. Check HELP PROBE');
            h.probe=in_probe;
        end
        function h=set.scan(h,in_scan)
            assert(isa(in_scan,'huff.scan'), 'The input scan is not a SCAN class. Check HELP SCAN');
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
             assert(isa(in_window,'huff.window'),'The window input should be a WINDOW class. Check help WINDOW');
             h.window=in_window;
        end
    end
    
    %% get method
    methods
        function value=get.data(h)
           
            % checking we have all we need
            %assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(numel(h.scan)>0,'The SCAN parameter is not set.');
            assert(numel(h.apex)>0,'The APEX parameter is not set.');
            %assert(all(h.scan.z>0),'Cannot compute apodization for points behind the XY plane');
            
            % NONE APODIZATION
            if(h.window==huff.window.none)
                value=ones(h.scan.N_pixels,h.N_elements);
                return;
            end
            
            % STA APODIZATION (just the element closest to apex)
            if (h.window==huff.window.sta)
                assert(numel(h.probe)>0,'The PROBE parameter is not set.');
                dist=sqrt((h.probe.x-h.apex.x).^2+(h.probe.y-h.apex.y).^2+(h.probe.z-h.apex.z).^2);
                value=ones(h.scan.N_pixels,1)*double(dist==min(dist(:)));
                return;
            end
            
            % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
            %x_dist=abs(h.element_position_matrix(:,:,1) - h.origin(:,1)*ones(1,h.probe.N_elements));
            x_dist=abs(bsxfun(@minus,h.element_position_matrix(:,:,1),h.origin(:,1)));
            %y_dist=abs(h.element_position_matrix(:,:,2) - h.origin(:,2)*ones(1,h.probe.N_elements));    
            y_dist=abs(bsxfun(@minus,h.element_position_matrix(:,:,2),h.origin(:,2)));
            
            % computing aperture
            Aperture_x=(abs(h.scan.z)./h.f_number(1))*ones(1,h.N_elements);
            Aperture_y=(abs(h.scan.z)./h.f_number(2))*ones(1,h.N_elements);

            % SWITCH
            switch(h.window)
                % BOXCAR/FLAT/RECTANGULAR
                case huff.window.boxcar
                    value=double(x_dist<=Aperture_x/2).*...
                          double(y_dist<=Aperture_y/2);
                % HANNING
                case huff.window.hanning
                    value=double(x_dist<=Aperture_x/2).*(0.5 + 0.5*cos(2*pi*x_dist./Aperture_x)).*... 
                          double(y_dist<=Aperture_y/2).*(0.5 + 0.5*cos(2*pi*y_dist./Aperture_y));
                % HAMMING
                case huff.window.hamming
                    value=double(x_dist<=Aperture_x/2).*(0.53836 + 0.46164*cos(2*pi*x_dist./Aperture_x)).*...
                          double(y_dist<=Aperture_y/2).*(0.53836 + 0.46164*cos(2*pi*y_dist./Aperture_y));
                % TUKEY25        
                case huff.window.tukey25
                    roll=0.25;
                    value=(x_dist<=(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<=(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY50        
                case huff.window.tukey50
                    roll=0.50;
                    value=(x_dist<=(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<=(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY75        
                case huff.window.tukey75
                    roll=0.75;
                    value=(x_dist<=(Aperture_x/2*(1-roll))) + (x_dist>(Aperture_x/2*(1-roll))).*(x_dist<(Aperture_x/2)).*0.5.*(1+cos(2*pi/roll*(x_dist./Aperture_x-roll/2-1/2))).*...                               
                          (y_dist<=(Aperture_y/2*(1-roll))) + (y_dist>(Aperture_y/2*(1-roll))).*(y_dist<(Aperture_y/2)).*0.5.*(1+cos(2*pi/roll*(y_dist./Aperture_y-roll/2-1/2)));
                % TUKEY80        
                case huff.window.tukey80
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
        
        function value=get.origin(h)
            % checking we have all we need
            assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(numel(h.scan)>0,'The SCAN parameter is not set.');
            assert(numel(h.apex)>0,'The APEX parameter is not set.');
            
            % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
            if isinf(h.apex.distance)
                origin(:,1)=h.scan.x+h.scan.z*tan(h.apex.azimuth)+h.scan.z*tan(h.tilt(1));
                origin(:,2)=h.scan.y+h.scan.z*tan(h.apex.elevation)+h.scan.z*tan(h.tilt(2));
            else
                origin(:,1)=h.apex.x-h.apex.z*(h.scan.x-h.apex.x)./(h.scan.z-h.apex.z)+h.scan.z*tan(h.tilt(1));
                origin(:,2)=(h.scan.y-h.apex.y)./(h.scan.x-h.apex.x).*(origin(:,1)-h.apex.x)+h.scan.y+h.scan.z*tan(h.tilt(2));
            end
            origin(:,3)=0;
            origin(isnan(origin))=0; % solve divisions by 0
            
            value=origin;
        end
        
        function value=get.propagation_distance(h)
            value=sqrt(sum(([h.scan.x-h.origin(:,1) h.scan.y-h.origin(:,2) h.scan.z-h.origin(:,3)]).^2,2));
        end
        
        function value=get.element_position_matrix(h)
            if isempty(h.sequence)
                assert(numel(h.probe)>0,'The PROBE parameter is not set.');
                value(:,:,1)=ones(h.scan.N_pixels,1)*(h.probe.x.'); 
                value(:,:,2)=ones(h.scan.N_pixels,1)*(h.probe.y.'); 
                value(:,:,3)=ones(h.scan.N_pixels,1)*(h.probe.z.'); 
            else
                % the transmit aperture in a synthetic scenario (CPWC, DW)
                % is dependant on the pixel location according to 10.1109/TUFFC.2015.007183 
                % (assuming flat apertures, not accurate for curvilinear probes)
                for n=1:length(h.sequence)
                    if isinf(h.sequence(n).source.distance)
                        element(:,1)=h.scan.x+h.scan.z*tan(h.sequence(n).source.azimuth);
                        element(:,2)=h.scan.y+h.scan.z*tan(h.sequence(n).source.elevation);
                    else
                        element(:,1)=h.sequence(n).source.x-h.sequence(n).source.z*(h.scan.x-h.sequence(n).source.x)./(h.scan.z-h.sequence(n).source.z);
                        element(:,2)=(h.scan.y-h.sequence(n).source.y)./(h.scan.x-h.sequence(n).source.x).*(element(:,1)-h.sequence(n).source.x)+h.scan.y;
                    end
                    element(:,3)=0;
                    element(isnan(element))=0; % solve divisions by 0
                    
                    value(:,n,1)=element(:,1);
                    value(:,n,2)=element(:,2);
                    value(:,n,3)=element(:,3);
                end
            end
        end
        
        function value=get.N_elements(h)
            if isempty(h.sequence)
                assert(numel(h.probe)>0,'The PROBE parameter is not set.');
                value=h.probe.N_elements; 
            else
                value=length(h.sequence);
            end
        
        end
    end
end