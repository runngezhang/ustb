classdef stack < process
%STACK   Matlab implementation of scanline stacking
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/05/01$

    properties
        name='Stack scanlines MATLAB';  % process name
        version='v1.0.0';                    % version
    end
    
    methods
        function out_dataset = go(h)
            % check input
            assert(isa(h.beamformed_data(1).scan,'uff.linear_scan')||isa(h.beamformed_data(1).scan,'uff.sector_scan'),'Stack only works with LINEAR_SCAN and SECTOR_SCAN');
            [Nrx Ntx]=size(h.beamformed_data);
            N=Nrx*Ntx;
            tools.workbar();
            switch class(h.beamformed_data(1).scan)
                case 'uff.linear_scan'
                    for nrx=1:Nrx
                        out_dataset(nrx,1)=uff.beamformed_data();
                        x_axis=h.beamformed_data(nrx,1).scan.x(1);
                        z_axis=h.beamformed_data(nrx,1).scan.z;
                        for ntx=1:Ntx
                            n=(nrx-1)*Ntx+ntx;
                            tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                            x_axis(ntx)=h.beamformed_data(nrx,ntx).scan.x(1);
                            out_dataset(nrx,1).data=[out_dataset(nrx,1).data; h.beamformed_data(nrx,ntx).data];
                        end
                        out_dataset(nrx,1).scan=uff.linear_scan(x_axis.',z_axis);
                    end
                case 'uff.sector_scan'
                    for nrx=1:Nrx
                        out_dataset(nrx,1)=uff.beamformed_data();
                        azimuth_axis=h.beamformed_data(nrx,1).scan.azimuth_axis(1);
                        depth_axis=h.beamformed_data(nrx,1).scan.depth_axis;
                        for ntx=1:Ntx
                            n=(nrx-1)*Ntx+ntx;
                            tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                            azimuth_axis(ntx)=h.beamformed_data(nrx,ntx).scan.azimuth_axis(1);
                            out_dataset(nrx,1).data=[out_dataset(nrx,1).data; h.beamformed_data(nrx,ntx).data];
                        end
                        out_dataset(nrx,1).scan=uff.sector_scan(azimuth_axis.',depth_axis);
                    end
            end
            tools.workbar(1);
        end
    end
end