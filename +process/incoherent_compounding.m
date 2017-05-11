classdef incoherent_compounding < process
    %INCOHERENT_COMPOUNDING   Matlab implementation of incoherent compounding
    %
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %            Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Last updated: 2017/05/11$
    
    %% constructor
    methods (Access = public)
        function h=incoherent_compounding()
            h.name='Incoherent compounding MATLAB';
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
                    out_data.data = abs(out_data.data);
                    % loop over beamformed dataset
                    N=Nrx*Ntx;
                    tools.workbar();
                    for n=2:N
                        if mod(n,round(N/100))==2
                            tools.workbar(n/N,sprintf('%s (%s) on "%s"',h.name,h.version,h.dimension),'USTB');
                        end
                        out_data.data=out_data.data+abs(h.beamformed_data(n).data);
                    end
                    tools.workbar(1);
                case dimension.transmit
                    tools.workbar();
                    for n_rx = 1:Nrx
                        tools.workbar(n_rx/Nrx,sprintf('%s (%s) on "%s"',h.name,h.version,h.dimension),'USTB');
                        % declare & copy beamformed dataset
                        out_data(n_rx,1)=uff.beamformed_data(h.beamformed_data(1));
                        out_data(n_rx,1).data = zeros(size(out_data(n_rx,1).data));
                        for n_tx = 1:Ntx
                            out_data(n_rx,1).data = out_data(n_rx,1).data + abs(h.beamformed_data(n_rx,n_tx).data);
                        end
                    end
                    tools.workbar(1);
                case dimension.receive
                    tools.workbar();
                    for n_tx = 1:Ntx
                        tools.workbar(n_tx/Ntx,sprintf('%s (%s) on "%s"',h.name,h.version,h.dimension),'USTB');
                        % declare & copy beamformed dataset
                        out_data(1,n_tx)=uff.beamformed_data(h.beamformed_data(1));
                        out_data(1,n_tx).data = zeros(size(out_data(1,n_tx).data));
                        for n_rx = 1:Nrx
                            out_data(1,n_tx).data = out_data(1,n_tx).data + abs(h.beamformed_data(n_rx,n_tx).data);
                        end
                    end
                    tools.workbar(1);
                otherwise
                    error('Unknown dimension mode; check HELP dimension');
            end
        end
    end
end
