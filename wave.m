classdef wave 
%wave   Wave definition
%
%   See also WAVE, SOURCE, PHANTOM, PROBE, PULSE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/22 $

    %% public properties
    properties  (SetAccess = public)
        delay            % delay [s]
        apodization      % apodization [unitless]
        source           % source class
    end
    
    %% dependent properties
    properties  (Dependent)   
        N_elements         % number of elements 
    end
    
    %% constructor
    methods (Access = public)
        function h=wave(source_in)
            %WAVE   Constructor of WAVE class
            %
            %   Syntax:
            %   h = wave(source)
            %           source      % source class
            %
            %   See also WAVE, SOURCE, PHANTOM, PROBE, PULSE
            
            if nargin>0
               h.source=source_in;
            end
        end
    end
    
    %% plot methods
    methods
        function plot(h)
          figure(1); 
          subplot(1,2,1);
          plot(1:numel(h.delay),h.delay*1e6); grid on; axis tight;
          xlabel('element');
          ylabel('delay [\mus]');
          set(gca,'fontsize',14);
          ylim([min([h.delay*1e6; -1]) max([1; h.delay*1e6])]);
          title('Delays');
          
          subplot(1,2,2);
          plot(1:numel(h.apodization),h.apodization); grid on; axis tight;
          xlabel('element');
          ylabel('apodization');
          set(gca,'fontsize',14);
          ylim([0 1.1]);
          title('Apodization');
          set(gcf,'Position',[261   605   859   343]);
          
          pause(0.01);
          
        end
    end
    
    %% set methods
    methods  
        function h=set.apodization(h,in_apodization)
            assert(size(in_apodization,2)==1, 'The apodization should be a column vector [N_elements 1]');
            h.apodization=in_apodization;
        end
        function h=set.delay(h,in_delay)
            assert(size(in_delay,2)==1, 'The delays should be a column vector [N_elements 1]');
            h.delay=in_delay;
        end
    end
    
    %% get methods
    methods  
        function value=get.N_elements(h)
            value=numel(h.delay);
        end
    end
    
end