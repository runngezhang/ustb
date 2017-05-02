classdef max < process
%MAX   Matlab implementation of max of beamformed values
%
%   authors: Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/05/02$

   
    %% constructor
    methods (Access = public)
        function h=max()
            h.name='Maximum value MATLAB'   
            h.reference= 'www.ustb.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};    
            h.version='v1.0.1';
        end
    end

    methods
        function out_data=go(h)
            % declare & copy beamformed dataset
            out_data=uff.beamformed_data(h.beamformed_data(1));
            
            % initialize
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
