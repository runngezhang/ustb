classdef coherence_factor < adaptive_beamformers.adaptive_beamformer
    %COHERENCE_FACTOR is a subclass of the adaptive_beamformer class
    %
    %   It implements the Mallart-Fink coherence factor from 
    %   R. Mallart and M. Fink, Adaptive focusing in scattering media through sound-speed inhomogeneities:
    %   The van Cittert Zernike approach and focusing criterion, J. Acoust. Soc. Am., vol. 96, no. 6, pp. 3721-3732, 1994
    %
    %   It also serves as an example on how to implement adaptive beamformers
    %   in the USTB.
    %
    %   authors: Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Date: 2017/04/07$
    
    % The coherence factor have no paramteres that the user needs to set,
    % thus the properties list can be left blank. For an example on a
    % method with parameters see the phase_coherence_factor
    properties
        
    end
    
    methods
        function CF_image = go(h,bmf)
            assert(strcmp(class(bmf),'beamformer'),'Input have to be a beamformer object');
            
            % First let's make the DAS image that the coherence factor will
            % be multiplied with.
            das = sum(h.data_cube,3);
            
            % A buffer to hold the calculated CF values
            % The buffer has the same size as the image
            cf = zeros(size(h.data_cube,1),size(h.data_cube,2));
            
            % Check if F number is zero for both transmit and receive
            % If so we can assume that the whole data_cube have nonzero
            % data so we don't have to index it.
            if h.fnumber_is_zero(bmf)
                
                % The CF is defined as the coherent sum divided
                % by the incoherent sum across the aperture.
                channels    = bmf.channel_data.N_elements;
                coherent    = abs(sum(h.data_cube,3)).^2;
                incoherent  = sum(abs(h.data_cube).^2,3);
                cf          = coherent./incoherent./channels;
                cf(isnan(cf)) = 0;
                CF_image = cf.*das;
            else
                % I'll fix this later
                % elseif strcmp(apod_type,'rx_apod')
                %     for zs = 1:size(h.data_cube,1)
                %         for xs = 1:size(h.data_cube,2)
                %             % Tx Rx apod
                %             valid_channels = h.data_cube(zs,xs,logical(squeeze(rx_apo(zs,xs,:))));
                %             coherent=abs(sum(valid_channels)).^2;
                %             incoherent=sum(abs(valid_channels.^2));
                %             cf(zs,xs) = coherent/incoherent/length(valid_channels);
                %         end
                %     end
                %     cf(cf<1*10^-6) = 10^-6;
                %     CF_image = das.*cf;
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
