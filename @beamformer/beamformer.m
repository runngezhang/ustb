classdef beamformer
%BEAMFORMER   beamformer definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Date: 2017/03/10$

    %% public properties
    properties  (SetAccess = public)
        channel_data         % CHANNEL_DATA class
        receive_apodization  % APODIZATION class
        transmit_apodization % APODIZATION class
        scan                 % collection of SCAN classes
    end
    
    %% optional properties
    properties  (SetAccess = public)
        pulse                % PULSE class
    end
    
    %% constructor
    methods (Access = public)
        function h=beamformer(bmf)
            %beamformer   Constructor of beamformer class
            %
            %   Syntax:
            %   h = beamformer()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE                      
            
            h.receive_apodization=uff.apodization();  % APODIZATION class
            h.transmit_apodization=uff.apodization(); % APODIZATION class
            
            % If called with another beamformer object as parameter
            % copy the properties
            if nargin > 0 && ~isempty(bmf)
                h.channel_data = bmf.channel_data;
                h.receive_apodization = bmf.receive_apodization;
                h.transmit_apodization = bmf.transmit_apodization;
                h.scan = bmf.scan;
            end
        end
    end
    
    methods 
        function out_dataset = go(h,postprocess)
            out_dataset = matlab_delay_base(h);
            
            if nargin == 1
                if numel(out_dataset) > 1
                    no_postprocess_warning()
                end
            else
                out_dataset = postprocess.go(out_dataset);
            end
        end
    end
%     %% set methods
%     methods  
%         function out_dataset=go(h,implementation,postprocess)
%             
%             % checking we have all we need
%             assert(numel(h.channel_data)>0,'The channel_data parameter is not set.');
%             assert(numel(h.scan)>0,'The SCAN parameter is not set.');
%             
%             %% beamforming
%             if ~exist('implementation')||isempty(implementation)
%                 inter_dataset=h.matlab();
%             elseif isa(implementation,'function_handle') % If implementation is a function handle, call it!
%                 inter_dataset=implementation();
%             elseif isobject(implementation)              % If it is a object, check if it is a adaptive beamformer 
%                 s = superclasses(implementation);        % subclass and call it!
%                 if strcmp(s{1,1},'adaptive_beamformers.adaptive_beamformer')
%                     inter_dataset=h.matlab(implementation);
%                 else
%                     error('Sorry! The adaptive beamformer have to be a subclass of the adaptive_beamformer class!');
%                 end
%             else
%                 error('Input for beamformer.go have to be a function_handle or a adaptive_beamformer object.');
%             end
% 
%             %% postprocess
%             if ~exist('postprocess')||isempty(postprocess)
%                 out_dataset=inter_dataset;
%             else
%                 out_dataset=postprocess(inter_dataset);
%             end
%             
%             %% If adaptive beamforming have been used, we want to save the parameters
%             if isobject(implementation)
%                 % Saving the adaptive beamformer so that we have the
%                 % parameters used
%                 out_dataset.adaptive_beamformer = implementation;
%             end
%         end
%     end
    
    %% set methods
    methods  
        function h=set.pulse(h,in_pulse)
            assert(isa(in_pulse,'uff.pulse'), 'The input is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.receive_apodization(h,in_apodization)
            assert(isa(in_apodization,'uff.apodization'), 'The input is not a APODIZATION class. Check HELP APODIZATION.');
            h.receive_apodization=in_apodization;
        end        
        function h=set.transmit_apodization(h,in_apodization)
            assert(isa(in_apodization,'uff.apodization'), 'The input is not a APODIZATION class. Check HELP APODIZATION.');
            h.transmit_apodization=in_apodization;
        end           
        function h=set.channel_data(h,in_channel_data)
            assert(isa(in_channel_data,'uff.channel_data'), 'The input is not a CHANNEL_DATA class. Check HELP CHANNEL_DATA.');
            h.channel_data=in_channel_data;
        end   
    end
    
    
end