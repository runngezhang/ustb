classdef apodization < uff
    %APODIZATION   UFF class to hold apodization data
    %   APODIZATION contains data to define transmit, receive & synthetic 
    %   beams. Different parameters are needed depending on the use.
    %
    %   Properties:
    %         probe               % UFF.PROBE class (needed for transmit & receive apodization)
    %         focus               % UFF.SCAN class (needed for transmit, receive & synthetic apodization)       
    %         sequence            % collection of UFF.WAVE classes (needed for synthetic apodizaton)
    %         
    %         window              % UFF.WINDOW class, default uff.window.noen
    %         f_number            % F-number [Fx Fy] [unitless unitless] 
    %         
    %         origo               % POINT class
    %         tilt                % tilt angle [azimuth elevation] [rad rad] 
    %
    %   Example:
    %         apo = uff.apodization();
    %
    %   See also UFF.CHANNEL_DATA, UFF.BEAMFORMED_DATA, UFF.SCAN
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %   $Last updated: 2017/06/09$

    %% public properties
    properties  (Access = public)
        probe                           % UFF.PROBE class (needed for transmit & receive apodization)
        focus                           % UFF.SCAN class (needed for transmit, receive & synthetic apodization)       
        sequence                        % collection of UFF.WAVE classes (needed for synthetic apodizaton)

        f_number  = [1 1]               % F-number [Fx Fy] [unitless unitless] 
        window    = uff.window.none     % UFF.WINDOW class, default uff.window.none
        
        origo     = uff.point('xyz',[0, 0, -Inf]);  % POINT class
        tilt      = [0 0]                           % tilt angle [azimuth elevation] [rad rad] 
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
        function h=apodization(varargin)
            h = h@uff(varargin{:});
        end
    end
    
    %% set methods
    methods  
        function h=set.origo(h,in_source)
            assert(isa(in_source,'uff.point'), 'The input origo is not a POINT class. Check HELP POINT');
            h.origo=in_source;
        end
        function h=set.probe(h,in_probe)
            assert(isa(in_probe,'uff.probe'), 'The input probe is not a PROBE class. Check HELP PROBE');
            h.probe=in_probe;
        end
        function h=set.focus(h,in_scan)
            assert(isa(in_scan,'uff.scan'), 'The input focus is not a SCAN class. Check HELP SCAN');
            h.focus=in_scan;
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
             assert(isa(in_window,'uff.window'),'The window input should be a WINDOW class. Check help WINDOW');
             h.window=in_window;
        end
    end
    
    methods
        function value=rectangular(h,dist,aperture)
            value=double(dist<=aperture/2);
        end
        function value=hanning(h,dist,aperture)
            value=double(dist<=aperture/2).*(0.5 + 0.5*cos(2*pi*dist./aperture));
        end
        function value=hamming(h,dist,aperture)
            value=double(dist<=aperture/2).*(0.53836 + 0.46164*cos(2*pi*dist./aperture));
        end
        function value=tukey(h,dist,aperture,roll)
            value=(dist<=(aperture/2*(1-roll))) + (dist>(aperture/2*(1-roll))).*(dist<(aperture/2)).*0.5.*(1+cos(2*pi/roll*(dist./aperture-roll/2-1/2)));
        end
    end
    
    %% get method
    methods
        function value=get.data(h)
           
            % checking we have all we need
            if isempty(h.focus)
                h.focus=uff.scan('xyz',[0 0 0]);
            end
            if isempty(h.origo)
                h.origo=uff.point('xyz',[0 0 0]);
            end
            
            % NONE APODIZATION
            if(h.window==uff.window.none)
                value=ones(h.focus.N_pixels,h.N_elements);
                return;
            end
            
            % STA APODIZATION (just the element closest to origo)
            if (h.window==uff.window.sta)
                assert(numel(h.probe)>0,'The PROBE parameter is not set.');
                dist=sqrt((h.probe.x-h.origo.x).^2+(h.probe.y-h.origo.y).^2+(h.probe.z-h.origo.z).^2);
                value=ones(h.focus.N_pixels,1)*double(dist==min(dist(:)));
                return;
            end
            
            % compute lateral distance 
            x_dist=abs(bsxfun(@minus,h.element_position_matrix(:,:,1),h.origin(:,1)));
            y_dist=abs(bsxfun(@minus,h.element_position_matrix(:,:,2),h.origin(:,2)));
            
            % computing aperture
            Aperture_x=(abs(h.focus.z)./h.f_number(1))*ones(1,h.N_elements);
            Aperture_y=(abs(h.focus.z)./h.f_number(2))*ones(1,h.N_elements);

            % SWITCH
            switch(h.window)
                % BOXCAR/FLAT/RECTANGULAR
                case uff.window.boxcar
                    value=h.rectangular(x_dist,Aperture_x).*h.rectangular(y_dist,Aperture_y);
                % HANNING
                case uff.window.hanning
                    value=h.hanning(x_dist,Aperture_x).*h.hanning(y_dist,Aperture_y);
                % HAMMING
                case uff.window.hamming
                    value=h.hamming(x_dist,Aperture_x).*h.hamming(y_dist,Aperture_y);
                % TUKEY25        
                case uff.window.tukey25
                    roll=0.25;
                    value=h.tukey(x_dist,Aperture_x,roll).*h.tukey(y_dist,Aperture_y,roll);
                % TUKEY50        
                case uff.window.tukey50
                    roll=0.50;
                    value=h.tukey(x_dist,Aperture_x,roll).*h.tukey(y_dist,Aperture_y,roll);
                % TUKEY75        
                case uff.window.tukey75
                    roll=0.75;
                    value=h.tukey(x_dist,Aperture_x,roll).*h.tukey(y_dist,Aperture_y,roll);
                % TUKEY80        
                case uff.window.tukey80
                    roll=0.80;
                    value=h.tukey(x_dist,Aperture_x,roll).*h.tukey(y_dist,Aperture_y,roll);
                otherwise
                    error('Unknown apodization type!');
            end
            
%             % edge smoothing
%             if(size(apo,2)>(2*(beam.smoothing+1)))
%                 % compute edge smoothing mask
%                 mask=0.5-0.5*cos((0:beam.smoothing)/beam.smoothing*pi);  % vector mask
%                 mask=(ones(h.focus.pixels,1)*mask(2:(beam.smoothing+1))); % matrix mask
%                 
%                 % modifying apodization edges
%                 apo(:,1:beam.smoothing)=apo(:,1:beam.smoothing).*mask;
%                 apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events)=apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events).*fliplr(mask);
%             end
        end
        
        function value=get.origin(h)
            % checking we have all we need
            %assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(h.focus.N_pixels>0,'The focus parameter is not set.');
            assert(numel(h.origo)>0,'The origo parameter is not set.');
            
            % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
            if isinf(h.origo.distance)
                origin(:,1)=h.focus.x+h.focus.z*tan(h.origo.azimuth)+h.focus.z*tan(h.tilt(1));
                origin(:,2)=h.focus.y+h.focus.z*tan(h.origo.elevation)+h.focus.z*tan(h.tilt(2));
            else
                origin(:,1)=h.origo.x-h.origo.z*(h.focus.x-h.origo.x)./(h.focus.z-h.origo.z)+h.focus.z*tan(h.tilt(1));
                origin(:,2)=(h.focus.y-h.origo.y)./(h.focus.x-h.origo.x).*(origin(:,1)-h.origo.x)+h.focus.y+h.focus.z*tan(h.tilt(2));
            end
            origin(:,3)=0;
            origin(isnan(origin))=0; % solve divisions by 0
            
            value=origin;
        end
        
        function value=get.propagation_distance(h)
            value=sqrt(sum(([h.focus.x-h.origin(:,1) h.focus.y-h.origin(:,2) h.focus.z-h.origin(:,3)]).^2,2));
        end
        
        function value=get.element_position_matrix(h)
            if isempty(h.sequence)
                assert(numel(h.probe)>0,'The PROBE parameter is not set.');
                value(:,:,1)=ones(h.focus.N_pixels,1)*(h.probe.x.'); 
                value(:,:,2)=ones(h.focus.N_pixels,1)*(h.probe.y.'); 
                value(:,:,3)=ones(h.focus.N_pixels,1)*(h.probe.z.'); 
            else
                % the transmit aperture in a synthetic scenario (CPWC, DW)
                % is dependant on the pixel location according to 10.1109/TUFFC.2015.007183 
                % (TODO: this is assuming a flat aperture, not accurate for curvilinear probes)
                for n=1:length(h.sequence)
                    if isinf(h.sequence(n).source.distance)
                        element(:,1)=h.focus.x+h.focus.z*tan(h.sequence(n).source.azimuth);
                        element(:,2)=h.focus.y+h.focus.z*tan(h.sequence(n).source.elevation);
                    else
                        element(:,1)=h.sequence(n).source.x-h.sequence(n).source.z*(h.focus.x-h.sequence(n).source.x)./(h.focus.z-h.sequence(n).source.z);
                        element(:,2)=h.sequence(n).source.y-h.sequence(n).source.z*(h.focus.y-h.sequence(n).source.y)./(h.focus.z-h.sequence(n).source.z);
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