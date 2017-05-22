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
            h.version='v1.0.3';
        end
    end
    
    properties (Access = public)
       dimension = dimension.both;          %Which "dimension" to sum over
    end

    methods
        function out_data=go(h)
            [N_pixels Nrx Ntx N_frames]=size(h.beamformed_data.data);
            
            switch h.dimension
                case dimension.both
                    out_data=uff.beamformed_data(h.beamformed_data); % ToDo: instead we should copy everything but the data
                    out_data.data=zeros(N_pixels, 1, 1, N_frames);
                    out_data.data=max(max(abs(h.beamformed_data.data),[],2),[],3);
                case dimension.transmit
                    out_data=uff.beamformed_data(h.beamformed_data); % ToDo: instead we should copy everything but the data
                    out_data.data=zeros(N_pixels, Nrx, 1, N_frames);
                    out_data.data=max(abs(h.beamformed_data.data),[],3);
                case dimension.receive
                    out_data=uff.beamformed_data(h.beamformed_data); % ToDo: instead we should copy everything but the data
                    out_data.data=zeros(N_pixels, 1, Ntx, N_frames);
                    out_data.data=max(abs(h.beamformed_data.data),[],2);
                otherwise
                    error('Unknown dimension mode; check HELP dimension');
            end
        end
    end 
end
