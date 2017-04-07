classdef phase_coherence_factor < adaptive_beamformers.adaptive_beamformer
    %PHASE_COHERENCE_FACTOR is a subclass of the adaptive_beamformer class
    %
    %    It implements the  Camacho-Fritsch Phase Coherence factor from
    %    http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4976281&tag=1
    %
    %   It also serves as an example on how to implement adaptive beamformers
    %   in the USTB.
    %
    %   authors: Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Date: 2017/04/07$
    
    %   This serves as an example on how to add parameters for the adaptive
    %   beamformer. These parameters could also have default values.
    properties
        gamma = 1 % In the article Comatcho-Fritz defines a parameter gamma,
                  % default value is 1.
    end
    
    methods
        function PCF_image = go(h,bmf)
            assert(strcmp(class(bmf),'beamformer'),'Input have to be a beamformer object');
            
            % First let's make the DAS image that the coherence factor will
            % be multiplied with.
            das = sum(h.data_cube,3);
            
            % Buffer for the "sf" values
            sf = zeros(size(h.data_cube,1),size(h.data_cube,2));
            
            % Check if F-number is zero for both transmit and receive
            % If so we can assume that the whole data_cube have nonzero
            % data so we don't have to index it.
            if h.fnumber_is_zero(bmf)
                
                % Please see the article referred for details on the
                % implementation
                sigma_0=pi/sqrt(3)
                phase=angle(h.data_cube);
                
                % auxiliary phase
                mask=phase<=0;
                A_phase=phase-pi;
                A_phase(mask)=phase(mask)+pi;
                
                % std deviation
                s_phase=std(phase,[],3);
                s_A_phase=std(A_phase,[],3);
                sf=bsxfun(@min,s_phase,s_A_phase);
                
                % APCF
                pcf=1-(h.gamma/sigma_0)*sf;
                pcf(pcf < 0) = 0;
                %mask=(pcf<0);
                pcf(isnan(pcf)) = 0;
                PCF_image=pcf.*das;
            else
                %                 for zs = 1:size(data_cube,1)
                %                     for xs = 1:size(data_cube,2)
                %                         phase=angle(data_cube(zs,xs,logical(squeeze(h.rx_apo(zs,xs,:)))));
                %                         % auxiliary phase
                %                         mask=phase<=0;
                %                         A_phase=phase-pi;
                %                         A_phase(mask)=phase(mask)+pi;
                %                         % std deviation
                %                         s_phase=std(phase,[],3);
                %                         s_A_phase=std(A_phase,[],3);
                %                         sf(zs,xs)=bsxfun(@min,s_phase,s_A_phase);
                %                     end
                %                 end
                error('CF only defined for f_number == 0 at the moment')
            end
        end
    end
    
    
    methods (Access = private)
        function result = fnumber_is_zero(h,bmf)
            % Check if F number is zero for both transmit and receive
            result = (sum(bmf.receive_apodization.f_number == 0) && ...
                ~isempty(bmf.transmit_apodization)) || ...
                sum(bmf.receive_apodization.f_number == 0) && ...
                sum(bmf.transmit_apodization.f_number == 0);
        end
    end
end