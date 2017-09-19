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
%   $Last updated: 2017/09/12$

    %% constructor
    methods (Access = public)
        function h=coherence_factor()
            h.name='Coherence Factor MATLAB';   
            h.reference= 'R. Mallart and M. Fink, Adaptive focusing in scattering media through sound-speed inhomogeneities: The van Cittert Zernike approach and focusing criterion, J. Acoust. Soc. Am., vol. 96, no. 6, pp. 3721-3732, 1994';                
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};    
            h.version='v1.0.5';
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
            
            % check dimensions
            if (h.dimension==dimension.receive) && (h.input.N_channels<2)
                error('Not enough channels to compute factor');
            end
            if (h.dimension==dimension.transmit) && (h.input.N_waves<2)
                error('Not enough waves to compute factor');
            end
            if (h.dimension==dimension.both) 
                if (h.input.N_channels<2)&&(h.input.N_waves>1)
                    warning('Not enough channels to compute factor. Changing dimension to dimension.transmit');
                    h.dimension = dimension.transmit;
                elseif (h.input.N_waves<2)&&(h.input.N_channels>1)
                    warning('Not enough waves to compute factor. Changing dimension to dimension.receive');
                    h.dimension = dimension.receive;
                elseif (h.input.N_waves<2)&&(h.input.N_channels<2)
                    error('Not enough waves and channels to compute factor');
                end
            end
                        
            % check if we have information about the receive apodization
            if (h.dimension~=dimension.transmit)
                if isempty(h.receive_apodization)||(h.receive_apodization.window==uff.window.none)
                    rx_apodization=ones(h.input.N_pixels,h.input.N_channels);
                else
                    h.receive_apodization.focus = h.input.scan;
                    rx_apodization=h.receive_apodization.data;
                end
            end
            
            % check if we have information about the transmit apodization
            if (h.dimension~=dimension.receive)
                if isempty(h.transmit_apodization)||(h.transmit_apodization.window==uff.window.none)||isempty(h.input.sequence)
                    tx_apodization=ones(h.input.N_pixels, 1, h.input.N_waves);
                else
                    h.transmit_apodization.focus = h.input.scan;
                    % calculate transmit apodization according to 10.1109/TUFFC.2015.007183
                    h.transmit_apodization.sequence=h.input.sequence;
                    tx_apodization=reshape(h.transmit_apodization.data,[h.input.N_pixels, 1, h.input.N_waves]);
                end
            end

            % declare output structure
            h.output=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data
            h.CF=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data

            switch h.dimension
                case dimension.both
                    active_elements=double(bsxfun(@times,tx_apodization,rx_apodization)>h.active_element_criterium);
                    coherent_sum=sum(sum(h.input.data,2),3);
                    incoherent_2_sum=sum(sum(abs(h.input.data).^2,2),3);
                    M=sum(sum(active_elements,2),3);
                case dimension.transmit
                    active_elements=double(tx_apodization>h.active_element_criterium);
                    coherent_sum=sum(h.input.data,3);
                    incoherent_2_sum=sum(abs(h.input.data).^2,3);
                    M=sum(active_elements,3);
                case dimension.receive
                    active_elements=double(rx_apodization>h.active_element_criterium);
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
    
    %% set methods
    methods
        function h=set.receive_apodization(h,in_apodization)
            assert(isa(in_apodization,'uff.apodization'), 'The input is not a UFF.APODIZATION class. Check HELP UFF.APODIZATION.');
            h.receive_apodization=in_apodization;
        end
        
        function h=set.transmit_apodization(h,in_apodization)
            assert(isa(in_apodization,'uff.apodization'), 'The input is not a UFF.APODIZATION class. Check HELP UFF.APODIZATION.');
            h.transmit_apodization=in_apodization;
        end
    end

end


