classdef coherence_factor < beamformer
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
        data_cube
        apo
    end
    
    
    
    methods
        %Constructor
        function h = coherence_factor(bmf)
            if nargin < 1
                bmf = [];
            end
            h = h@beamformer(bmf);
        end
        
        
        function out_dataset = go(h,postprocess)
            out_dataset = matlab_delay_base(h,@coherence_factor_implementation);
            
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
        function CF_image = coherence_factor_implementation(h)
            %assert(strcmp(class(bmf),'beamformer'),'Input have to be a beamformer object');
            
            % First let's make the DAS image that the coherence factor will
            % be multiplied with.
            das = sum(h.data_cube,3);
            
            % A buffer to hold the calculated CF values
            % The buffer has the same size as the image
            cf = zeros(size(h.data_cube,1),size(h.data_cube,2));
            
            % Check there is apodization on transmit or receive
            % If not we can assume that the whole data_cube have nonzero
            % data so we don't have to index it.
            if h.transmit_apodization.window == 'none' && h.receive_apodization.window == 'none'
                
                % The CF is defined as the coherent sum divided
                % by the incoherent sum across the aperture.
                channels    = size(h.data_cube,3);
                coherent    = abs(sum(h.data_cube,3)).^2;
                incoherent  = sum(abs(h.data_cube).^2,3);
                cf          = coherent./incoherent./channels;
                cf(isnan(cf)) = 0;
                CF_image = cf.*das;
            else %If there is apodization we have to run through the channel data,
                %% is there a better way???
                assert(h.transmit_apodization.window == 'boxcar' || h.transmit_apodization.window == 'none','Please use boxcar apodization for transmit when using coherence factor');
                assert(h.receive_apodization.window == 'boxcar','Please use boxcar apodization receive when using coherence factor');
                %apod_matrix = reshape(h.apo,h.scan.N_z_axis,h.scan.N_x_axis,h.channel_data.probe.N_elements);
                apod_matrix = reshape(h.apo,size(h.data_cube,1),size(h.data_cube,2),size(h.data_cube,3));
          
                for zs = 1:size(h.data_cube,1)
                    for xs = 1:size(h.data_cube,2)
                        valid_channels = h.data_cube(zs,xs,logical(squeeze(apod_matrix(zs,xs,:))));
                        coherent=abs(sum(valid_channels)).^2;
                        incoherent=sum(abs(valid_channels.^2));
                        if incoherent ~= 0 % We don't want NaN's
                            cf(zs,xs) = coherent/incoherent/length(valid_channels);
                        end
                    end
                end
                %cf(cf<1*10^-6) = 10^-6;
                CF_image = das.*cf;
                
            end
        end
    end
end
