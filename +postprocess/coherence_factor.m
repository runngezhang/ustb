classdef coherence_factor < postprocess
%COHERENCE_FACTOR   Matlab implementation of Mallart-Fink Coherence Factor
%
%   MATLAB implementation of Mallart-Fink coherence factor beamforming
%   method as described in the paper:
%
%   R. Mallart and M. Fink, "Adaptive focusing in scattering media through 
%   sound-speed inhomogeneities: The van Cittert Zernike approach and focusing 
%   criterion", J. Acoust. Soc. Am., vol. 96, no. 6, pp. 3721-3732, 1994
%
%   The implementation computes coherence either on transmit, receive, or
%   both.
%
%   implementers: Ole Marius Hoel Rindal <olemarius@olemarius.net>
%                 Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%
%   $Last updated: 2017/09/10$

    %% constructor
    methods (Access = public)
        function h=coherence_factor()
            h.name='Coherence Factor MATLAB';   
            h.reference= 'R. Mallart and M. Fink, Adaptive focusing in scattering media through sound-speed inhomogeneities: The van Cittert Zernike approach and focusing criterion, J. Acoust. Soc. Am., vol. 96, no. 6, pp. 3721-3732, 1994';                
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};    
            h.version='v1.0.4';
        end
    end
        
    %% Additional properties
    properties
        CF                                            % BEAMFORMED_DATA class with the computed coherent factor
        receive_apodization                           % APODIZATION class
        transmit_apodization                          % APODIZATION class
        active_element_criterium=0.16;                % value to decide whether an element is used or not. This value depends on the SNR so it must be adjusted on a case-by-case basis.
        dimension = dimension.both;                   % dimension class that specifies whether the process will run only on transmit, receive, or both.
    end

    methods
        function output = go(h)
            % check if we can skip calculation
            if h.check_hash()
                output= h.output; 
                return;
            end
           
            % check if we have information about apodization
            rx_apodization=ones(h.input.N_pixels,h.input.N_channels);
            tx_apodization=ones(h.input.N_pixels,h.input.N_waves);
            if ~isempty(h.transmit_apodization)&~isempty(h.receive_apodization)&~isempty(h.channel_data.probe)
                % receive
                if h.input.N_channels>1
                    h.receive_apodization.focus=h.input.scan;
                    h.receive_apodization.probe=h.channel_data.probe;
                    rx_apodization=h.receive_apodization.data();                
                end
                
                % transmit
                if h.input.N_waves>1
                    h.transmit_apodization.sequence = h.channel_data.sequence;
                    h.transmit_apodization.focus=h.input.scan;
                    h.transmit_apodization.probe=h.channel_data.probe;
                    tx_apodization=h.transmit_apodization.data();
                end
            else
                warning('Missing probe and apodization data; full aperture is assumed.');
            end
            tx_apodization=reshape(tx_apodization,[h.input.N_pixels, 1, h.input.N_waves]);
            apodization_matrix=bsxfun(@times,tx_apodization,rx_apodization);
            active_elements=double(apodization_matrix>h.active_element_criterium);
            
            % declare output structure
            h.output=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data
            h.CF=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data

            switch h.dimension
                case dimension.both
                    coherent_sum=sum(sum(h.input.data,2),3);
                    incoherent_2_sum=sum(sum(abs(h.input.data).^2,2),3);
                    M=sum(sum(active_elements,2),3);
                case dimension.transmit
                    coherent_sum=sum(h.input.data,3);
                    incoherent_2_sum=sum(abs(h.input.data).^2,3);
                    M=sum(active_elements,3);
                case dimension.receive
                    coherent_sum=sum(h.input.data,2);
                    incoherent_2_sum=sum(abs(h.input.data).^2,2);
                    M=sum(active_elements,2);
                otherwise
                    error('Unknown dimension mode; check HELP dimension');
            end
            % Coherent Factor
            h.CF.data = bsxfun(@times,abs(coherent_sum).^2./incoherent_2_sum,1./M); 
            % coherent factor image            
            h.output.data = h.CF.data .* coherent_sum;
            
            % pass reference
            output = h.output;
            
            % update hash
            h.save_hash();
        end   
    end
end



