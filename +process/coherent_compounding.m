classdef coherent_compounding < process
%COHERENT_COMPOUNDING   Matlab implementation of Coherent compounding
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/04/24$

    properties
        name='Coherent compounding MATLAB';  % process name
        version='v1.0.0';                    % version
    end

    methods
        function out_data=go(h)
            
            out_data=h.beamformed_data(1);
            tools.workbar();
            for n=2:length(h.beamformed_data(:))
                tools.workbar(n/length(h.beamformed_data(:)),sprintf('%s (%s)',h.name,h.version),'USTB');
                out_data.data=out_data.data+h.beamformed_data(n).data;
            end
            
        end
    end 
end
