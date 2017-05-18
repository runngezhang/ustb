classdef coherent_compounding < process
%COHERENT_COMPOUNDING   Matlab implementation of Coherent compounding
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/05/11$

    %% constructor
    methods (Access = public)
        function h=coherent_compounding()
            h.name='Coherent compounding MATLAB';   
            h.reference='www.ntnu.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};   
            h.version='v1.0.2';
        end
    end
    
    properties (Access = public)
       dimension = dimension.both;          %Which "dimension" to sum over
    end
    
    methods
        function out_data=go(h)
            [Nrx Ntx]=size(h.beamformed_data);
            
            switch h.dimension
                case dimension.both
                    % declare & copy beamformed dataset
                    out_data=uff.beamformed_data(h.beamformed_data(1));
                    % loop over beamformed dataset
                    N=Nrx*Ntx;
                    tools.workbar();
                    for n=2:N
                        if mod(n,round(N/100))==2
                            tools.workbar(n/N,sprintf('%s (%s) on "%s"',h.name,h.version,char(h.dimension)),'USTB');
                        end
                        out_data.data=out_data.data+h.beamformed_data(n).data;
                    end
                    tools.workbar(1);
                case dimension.transmit
                    tools.workbar();
                    for n_rx = 1:Nrx
                        tools.workbar(n_rx/Nrx,sprintf('%s (%s) on "%s"',h.name,h.version,char(h.dimension)),'USTB');
                        % declare & copy beamformed dataset
                        out_data(n_rx,1)=uff.beamformed_data(h.beamformed_data(1));
                        out_data(n_rx,1).data = zeros(size(out_data(n_rx,1).data));
                        for n_tx = 1:Ntx
                            out_data(n_rx,1).data = out_data(n_rx,1).data + h.beamformed_data(n_rx,n_tx).data;
                        end
                    end
                    tools.workbar(1);
                case dimension.receive
                    tools.workbar();
                    for n_tx = 1:Ntx
                        tools.workbar(n_tx/Ntx,sprintf('%s (%s) on "%s"',h.name,h.version,char(h.dimension)),'USTB');
                        % declare & copy beamformed dataset
                        out_data(1,n_tx)=uff.beamformed_data(h.beamformed_data(1));
                        out_data(1,n_tx).data = zeros(size(out_data(1,n_tx).data));
                        for n_rx = 1:Nrx
                            out_data(1,n_tx).data = out_data(1,n_tx).data + h.beamformed_data(n_rx,n_tx).data;
                        end
                    end
                    tools.workbar(1);
                otherwise
                    error('Unknown dimension mode; check HELP dimension');
            end
        end
    end
end
