classdef coherent_compounding < process
%COHERENT_COMPOUNDING   Matlab implementation of Coherent compounding
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/05/02$

    %% constructor
    methods (Access = public)
        function h=coherent_compounding()
            h.name='Coherent compounding MATLAB';   
            h.reference='www.ntnu.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};   
            h.version='v1.0.1';
        end
    end

    methods
        function out_data=go(h)
            % declare & copy beamformed dataset
            out_data=uff.beamformed_data(h.beamformed_data(1));
            
            % loop over beamformed dataset
            N=length(h.beamformed_data(:));
            tools.workbar();
            for n=2:N
                if mod(n,round(N/100))==2
                    tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                end
                out_data.data=out_data.data+h.beamformed_data(n).data;
            end
            tools.workbar(1);
        end
    end 
end
