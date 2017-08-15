classdef process < handle
    %PROCESS   process part of the beamforming chain
    %
    %   See also BEAMFORMER, CHANNEL_DATA, BEAMFORMED_DATA
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %            Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Date: 2017/04/28$
    
    %% public properties
    properties  (SetAccess = public)
        channel_data         % CHANNEL_DATA class
        receive_apodization  % APODIZATION class
        transmit_apodization % APODIZATION class
        scan                 % collection of SCAN classes
        beamformed_data      % collection of BEAMFORMED_DATA classes
    end
    
    %% Logistics
    properties  (SetAccess = public)
        name=''              % name of the process
        reference=''         % reference to the publication where it is disclossed
        implemented_by=''    % contact of the implementer/s
        version=''           % verion
    end
    
    %% constructor
    methods (Access = public)
        function h=process()
            %process   Constructor of process class
            %
            %   Syntax:
            %   h = process()
            %
            %   See also BEAMFORMER, CHANNEL_DATA, BEAMFORMED_DATA
            
            %h.channel_data=uff.channel_data();        % CHANNEL_DATA
            %h.receive_apodization=uff.apodization();  % APODIZATION class
            %h.transmit_apodization=uff.apodization(); % APODIZATION class
            %h.scan=uff.scan();                        % SCAN class
            %h.beamformed_data=uff.beamformed_data();  % BEAMFORMED_DATA class
            
            %             % If called with a beamformer object as parameter copy the properties
            %             if nargin > 0 && ~isempty(bmf) && isa(bmf,'beamformer')
            %                 h.channel_data = bmf.channel_data;
            %                 h.receive_apodization = bmf.receive_apodization;
            %                 h.transmit_apodization = bmf.transmit_apodization;
            %                 h.scan = bmf.scan;
            %             end
            %             % If called with a beamformed data object we copy the data
            %             if nargin > 1 && ~isempty(in_data) && isa(in_data,'uff.beamformed_data')
            %                 h.beamformed_data= in_data;
            %             end
        end
    end
    
    %% go method
    methods
        function out_data = go(h)
            % beamformed_data is a handle class. The equal operation (a=b)
            % passes the handle not a copy. To pass a copy the following
            % must be done
            
            [Nrx Ntx]=size(h.beamformed_data);
            for nrx=1:Nrx
                for ntx=1:Ntx
                    out_data(nrx,ntx) = uff.beamformed_data();
                    out_data(nrx,ntx).copy(h.beamformed_data(nrx,ntx));
                end
            end
            
            % handle classes allow us to limit memory use while updating
            % parameters within the class methods
        end
    end
    
    %% set methods
    methods
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
        function h=set.scan(h,in_scan)
            assert(isa(in_scan,'uff.scan'), 'The input is not a SCAN class. Check HELP SCAN.');
            h.scan=in_scan;
        end
        function h=set.beamformed_data(h,in_data)
            assert(isa(in_data,'uff.beamformed_data'), 'The input is not a BEAMFORMED_DATA class. Check HELP BEAMFORMED_DATA.');
            h.beamformed_data=in_data;
        end
    end
    
    %% Display methods
    methods
        function print_name(h)
            idx = 1;
            while idx+50 < length(h.name)
                if idx < 50
                    fprintf('Name:\t\t %s \n',h.name(idx:idx+49));
                else
                    fprintf('\t\t %s \n',h.name(idx:idx+49));
                end
                idx = idx+50;
            end
            if idx < 50
                fprintf('Name:\t\t %s \n',h.name(idx:end));
            else
                fprintf('\t\t %s \n',h.name(idx:end));
            end
        end
        
        function print_reference(h)
            idx = 1;
            while idx+50 < length(h.reference)
                if idx < 50
                    fprintf('Reference:\t %s \n',h.reference(idx:idx+49));
                else
                    fprintf('\t\t %s \n',h.reference(idx:idx+49));
                end
                idx = idx+50;
            end
            if idx < 50
                fprintf('Name:\t %s \n',h.reference(idx:end));
            else
                fprintf('\t\t %s \n',h.reference(idx:end));
            end
        end
        
        function print_implemented_by(h)
            idx = 1;
            while idx+50 < length(h.implemented_by)
                if idx < 50
                    fprintf('Implemented by:\t %s \n',h.implemented_by(idx:idx+49));
                else
                    fprintf('\t\t %s \n',h.implemented_by(idx:idx+49));
                end
                idx = idx+50;
            end
            if idx < 50
                fprintf('Name:\t %s \n',h.implemented_by(idx:end));
            else
                fprintf('\t\t %s \n',h.implemented_by(idx:end));
            end
        end
    end
    
end

