classdef max < process
%MAX   Matlab implementation of max of beamformed values
%
%   authors: Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/04/01$

    properties
        name='Maximum value MATLAB';  % process name
        version='v1.0.0';             % version
    end

    methods
        function out_data=go(h)
            out_data=h.beamformed_data(1);         
            out_data.data=abs(h.beamformed_data(1).data);
            N=length(h.beamformed_data(:));
            tools.workbar();
            for n=2:N
                if mod(n,round(N/100))==2
                    tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                end
                out_data.data=max(out_data.data, abs(h.beamformed_data(n).data));
            end
        end
    end 
end
