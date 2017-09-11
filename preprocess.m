classdef preprocess < process
    %PREPROCESS   preprocess part of the prcessing pipeline. Takes a
    % uff.channel_data structure and return another uff.channel_data
    %
    %   See also PROCESS, CHANNEL_DATA, BEAMFORMED_DATA
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %            Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Date: 2017/09/10$
    
    %% public properties
    properties  (Access = public)
        input                % CHANNEL_DATA class
        output               % CHANNEL_DATA class
    end

    %% Dependant properties
    properties  (Dependent)
        sampling_frequency   % sampling frequency [Hz]
    end

    
    %% set methods
    methods
        function h=set.input(h,in_channel_data)
            assert(isa(in_channel_data,'uff.channel_data'), 'The input is not a UFF.CHANNEL_DATA class. Check HELP UFF.CHANNEL_DATA.');
            h.input=in_channel_data;
        end
    end
    
    %% get methods
    methods
        function value=get.sampling_frequency(h)
            if isempty(h.input)
                error('The input channel data has not been asigned.');
            else
                value=h.input.sampling_frequency;
            end
        end
    end

end

