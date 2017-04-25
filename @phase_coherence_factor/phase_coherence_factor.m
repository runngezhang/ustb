classdef phase_coherence_factor < beamformer
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
        data_cube
        gamma = 1 % In the article Comatcho-Fritz defines a parameter gamma,
        % default value is 1.
        apo
    end
    
    methods
        %Constructor
        function h = phase_coherence_factor(bmf)
            if nargin < 1
                bmf = [];
            end
            h = h@beamformer(bmf);
        end
        
        function out_dataset = go(h,postprocess)
            out_dataset = matlab_delay_base(h,@phase_coherence_factor_implementation);
            
            if nargin == 1
                if numel(out_dataset) > 1
                    beamformer.no_postprocess_warning()
                end
            else
                out_dataset = postprocess.go(out_dataset);
            end
        end
    end
    
    methods
        function PCF_image = phase_coherence_factor_implementation(h)
            
            %  assert(strcmp(class(bmf),'beamformer'),'Input have to be a beamformer object');
            
            % First let's make the DAS image that the coherence factor will
            % be multiplied with.
            das = sum(h.data_cube,3);
            
            % Buffer for the "sf" values
            sf = zeros(size(h.data_cube,1),size(h.data_cube,2));
            
            % Check there is apodization on transmit or receive
            % If not we can assume that the whole data_cube have nonzero
            % data so we don't have to index it.
            if h.transmit_apodization.window == 'none' && h.receive_apodization.window == 'none'
                
                % Please see the article referred for details on the
                % implementation
                sigma_0=pi/sqrt(3);
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
                assert(h.transmit_apodization.window == 'boxcar' || h.transmit_apodization.window == 'none','Please use boxcar apodization for transmit when using coherence factor');
                assert(h.receive_apodization.window == 'boxcar','Please use boxcar apodization receive when using coherence factor');
                apod_matrix = reshape(h.apo,h.scan.N_z_axis,h.scan.N_x_axis,h.channel_data.probe.N_elements);
                
                sigma_0=pi/sqrt(3);
                for zs = 1:size(h.data_cube,1)
                    for xs = 1:size(h.data_cube,2)
                        phase=angle(h.data_cube(zs,xs,logical(squeeze(apod_matrix(zs,xs,:)))));
                        % auxiliary phase
                        mask=phase<=0;
                        A_phase=phase-pi;
                        A_phase(mask)=phase(mask)+pi;
                        % std deviation
                        s_phase=std(phase,[],3);
                        s_A_phase=std(A_phase,[],3);
                        sf(zs,xs)=bsxfun(@min,s_phase,s_A_phase);
                    end
                end
                % APCF
                pcf=1-(h.gamma/sigma_0)*sf;
                minimum_value_found = min(min(pcf(pcf > 0)));
                pcf(pcf < 0) = minimum_value_found;
                pcf(isnan(pcf)) = minimum_value_found;
                pcf(isinf(pcf)) = minimum_value_found;
                PCF_image=pcf.*das;
            end
        end
    end
end