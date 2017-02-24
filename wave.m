classdef wave 
%wave   Wave definition
%
%   See also WAVE, SOURCE, PHANTOM, PROBE, PULSE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/22 $

    %% public properties
    properties  (SetAccess = public)
        source           % source class
        probe            % probe class
        sound_speed      % reference speed of sound
    end
    
    %% dependent properties
    properties  (Dependent)   
        N_elements       % number of elements 
        delay            % delay [s]        
        apodization      % apodization [unitless]
    end
    
    %% constructor
    methods (Access = public)
        function h=wave(probe_in,source_in)
            %WAVE   Constructor of WAVE class
            %
            %   Syntax:
            %   h = wave(source)
            %           source      % source class
            %
            %   See also WAVE, SOURCE, PHANTOM, PROBE, PULSE
            
            if nargin>0
               h.probe=probe_in;
            end
            if nargin>1
               h.source=source_in;
            end            
        end
    end
    
    %% plot methods
    methods
        function fig_handle=plot(h)
          fig_handle=figure(); 

          % probe geometry
          x = [(h.probe.x-h.probe.width/2.*cos(h.probe.theta)).'; (h.probe.x+h.probe.width/2.*cos(h.probe.theta)).'; (h.probe.x+h.probe.width/2.*cos(h.probe.theta)).'; (h.probe.x-h.probe.width/2.*cos(h.probe.theta)).'];
          y = [(h.probe.y-h.probe.height/2.*cos(h.probe.phi)).'; (h.probe.y-h.probe.height/2.*cos(h.probe.phi)).'; (h.probe.y+h.probe.height/2.*cos(h.probe.phi)).'; (h.probe.y+h.probe.height/2.*cos(h.probe.phi)).'; ];
          z = zeros(size(x));
          c = linspace(0,1,h.probe.N_elements);
          
          %subplot(1,2,1);
          % draw flatten elements
          fill3(x*1e3,y*1e3,z*1e3,c); grid on; axis equal tight; hold on;
          % draw delays
          plot3(h.probe.x*1e3,h.probe.y*1e3,h.delay*1e6,'r.'); grid on; axis tight;
          xlabel('x [mm]');
          ylabel('y [mm]');
          zlabel('delay [\mus]');
          set(gca,'fontsize',14);
          title('Delays');
          
%           subplot(1,2,2);
%           % draw flatten elements
%           fill3(x*1e3,y*1e3,z*1e3,c); grid on; axis equal tight; hold on;
%           % draw apodization
%           plot3(h.probe.x*1e3,h.probe.y*1e3,h.apodization,'r.'); grid on; axis tight;
%           xlabel('x [mm]');
%           xlabel('y [mm]');
%           ylabel('Apodization');
%           set(gca,'fontsize',14);
%           title('Apodization');
          
        end
    end
    
    %% set methods
    methods  
%         function h=set.apodization(h,in_apodization)
%             assert(size(in_apodization,2)==1, 'The apodization should be a column vector [N_elements 1]');
%             assert(size(in_apodization,1)==h.N_elements, 'The apodization size does not match the probe definition');
%             h.apodization=in_apodization;
%         end
        function h=set.source(h,in_source)
            assert(strcmp(class(in_source),'source'), 'The source is not a SOURCE class. Check HELP SOURCE');
            h.source=in_source;
        end
        function h=set.probe(h,in_probe)
            assert(strcmp(class(in_probe),'probe'), 'The probe is not a PROBE class. Check HELP PROBE');
            h.probe=in_probe;
        end
    end
    
    %% get methods
    methods  
        function value=get.N_elements(h)
            value=h.probe.N_elements;
        end
        function value=get.delay(h)
            if ~isinf(h.source.distance)
                dst=sqrt((h.probe.x-h.source.x).^2+(h.probe.y-h.source.y).^2+(h.probe.z-h.source.z).^2);
                value=h.source.distance/h.sound_speed-dst/h.sound_speed;
            else
                value=h.probe.x*sin(h.source.azimuth)/h.sound_speed+h.probe.y*sin(h.source.elevation)/h.sound_speed;
            end
        end
        function value=get.apodization(h)
            % for the moment no apodization is programed
            value=ones(size(h.probe.x));
        end
    end
    
end