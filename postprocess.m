classdef postprocess < process
    %POSTPROCESS   postprocess part of the processing pipeline. Takes a
    % uff.beamformed_data structure and returns another uff.beamformed_data
    %
    %   See also PROCESS, CHANNEL_DATA, BEAMFORMED_DATA
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %            Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Date: 2017/09/10$
    
    %% public properties
    properties  (Access = public)
        input                % UFF.BEAMFORMED_DATA class
        output               % UFF.BEAMFORMED_DATA class
    end
    
    %% set methods
    methods
        function h=set.input(h,in_beamformed_data)
            assert(isa(in_beamformed_data,'uff.beamformed_data'), 'The input is not a UFF.BEAMFORMED_DATA class. Check HELP UFF.BEAMFORMED_DATA.');
            h.input=in_beamformed_data;
        end
    end
end

