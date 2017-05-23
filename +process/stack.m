classdef stack < process
%STACK   Matlab implementation of scanline stacking
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2017/05/02$

        
    %% constructor
    methods (Access = public)
        function h=stack()
            h.name='Stack scanlines MATLAB';          % name of the process
            h.reference='www.ustb.no';                % reference to the publication where it is disclossed
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};    
            h.version='v1.0.1';     
        end
    end
    
    methods
        function out_dataset = go(h)
            % check input
            assert(isa(h.beamformed_data(1).scan,'uff.linear_scan')||isa(h.beamformed_data(1).scan,'uff.sector_scan'),'Stack only works with LINEAR_SCAN and SECTOR_SCAN');
            [N_pixels Nrx Ntx N_frames]=size(h.beamformed_data.data);
            N=Nrx*Ntx;
            tools.workbar();
            switch class(h.beamformed_data(1).scan)
                case 'uff.linear_scan'
                    % declare
                    out_dataset=uff.beamformed_data(h.beamformed_data);     % ToDo: shouldn't copy the data
                    out_dataset.data=zeros(N_pixels*Ntx,Nrx,1,N_frames);
    
                    % loop over scans
                    depth_axis=h.beamformed_data.scan(1).z;
                    azimuth_axis=[];
                    for ntx=1:Ntx
                        tools.workbar(ntx/Ntx,sprintf('%s (%s)',h.name,h.version),'USTB');
                        azimuth_axis=[azimuth_axis h.beamformed_data.scan(ntx).x_axis];
                        out_dataset.data((1:N_pixels)+N_pixels*(ntx-1),:,1,:)=h.beamformed_data.data(:,:,ntx,:);
                    end
                    out_dataset.scan=uff.linear_scan(azimuth_axis.',depth_axis);
                case 'uff.sector_scan'
                    out_dataset=uff.beamformed_data(h.beamformed_data); % ToDo: shouldn't copy the data
                    out_dataset.data=zeros(N_pixels*Ntx,Nrx,1,N_frames);
    
                    % loop over scans
                    depth_axis=h.beamformed_data.scan(1).depth_axis;
                    azimuth_axis=[];
                    for ntx=1:Ntx
                        tools.workbar(ntx/Ntx,sprintf('%s (%s)',h.name,h.version),'USTB');
                        azimuth_axis=[azimuth_axis h.beamformed_data.scan(ntx).azimuth_axis];
                        out_dataset.data((1:N_pixels)+N_pixels*(ntx-1),:,1,:)=h.beamformed_data.data(:,:,ntx,:);
                    end
                    out_dataset.scan=uff.sector_scan(azimuth_axis.',depth_axis);
                    
            end
            tools.workbar(1);
        end
    end
end