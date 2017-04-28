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
    
    %% go method
    methods 
        function out_dataset = go(h,process_list)
            if nargin == 1
               proc=process.das_matlab();
               proc.channel_data=h.channel_data;
               proc.receive_apodization=h.receive_apodization;
               proc.transmit_apodization=h.transmit_apodization;
               proc.scan=h.scan;
               out_dataset = proc.go();
            else
               process_list{1}.channel_data=h.channel_data;
               process_list{1}.receive_apodization=h.receive_apodization;
               process_list{1}.transmit_apodization=h.transmit_apodization;
               process_list{1}.scan=h.scan;
               out_dataset = process_list{1}.go();
               for n=2:length(process_list)
                   process_list{n}.channel_data=h.channel_data;
                   process_list{n}.receive_apodization=h.receive_apodization;
                   process_list{n}.transmit_apodization=h.transmit_apodization;
                   process_list{n}.scan=h.scan;
                   process_list{n}.beamformed_data=out_dataset;
                   out_dataset = process_list{n}.go();
               end
            end
                
%             if nargin == 2 %Only one postprocess, this is assumed to be second
%                 postprocess_2 = postprocess_1;
%                 out_dataset = matlab_delay_interdataset(h,postprocess.delay_and_sum);
%             elseif nargin == 3
%                 out_dataset = matlab_delay_interdataset(h,postprocess_1);
%             end
%             
%             if nargin == 1
%                 if numel(out_dataset) > 1
%                     no_postprocess_warning()
%                 end
%             else
%                 postprocess_2.bmf = h;
%                 out_dataset = postprocess_2.go(out_dataset);
%             end
        end
    end
    
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