classdef phase_coherence_factor < postprocess
%PHASE_COHERENCE_FACTOR   Matlab implementation of Camacho-Fritsch Phase Coherence Factor
%
%   MATLAB implementation of Camacho-Fritsch Phase Coherence Factor beamforming
%   method as described in the paper:
%
%   J. Camacho and C. Fritsch, "Phase coherence imaging of grained materials," 
%   in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
%   vol. 58, no. 5, pp. 1006-1015, May 2011.
%
%   The implementation computes coherence either on transmit, receive, or
%   both.
%
%   implementers: Ole Marius Hoel Rindal <olemarius@olemarius.net>
%                 Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%
%   $Last updated: 2017/05/22$

        %% constructor
    methods (Access = public)
        function h=phase_coherence_factor()
            h.name='Phase Coherence Factor MATLAB';   
            h.reference='J. Camacho and C. Fritsch, Phase coherence imaging of grained materials, in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, vol. 58, no. 5, pp. 1006-1015, May 2011.';                
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};    
            h.version='v1.0.2';
        end
    end

    %% Additional properties
    properties
        gamma=1;                                      % mixing ratio
        sigma_0=pi/sqrt(3);                           % reference phase value
        receive_apodization                           % APODIZATION class
        transmit_apodization                          % APODIZATION class
        FCA                                           % BEAMFORMED_DATA class with the computed absolute phase coherence factor
        FCC                                           % BEAMFORMED_DATA class with the computed complex phase coherence factor
    end
    
    properties (Access = public)
       dimension = dimension.both;          %Which "dimension" to sum over
    end
    
    methods
        function output=go(h)
            % check if we can skip calculation
            if h.check_hash()
                output = h.output; 
                return;
            end            
           
            % check if we have information about apodization
            if isempty(h.receive_apodization)||(h.receive_apodization.window==uff.window.none)
                rx_apodization=ones(h.input.N_pixels,h.input.N_channels);
            else
                h.receive_apodization.focus = h.input.scan;
                rx_apodization=h.receive_apodization.data;
            end
            if isempty(h.transmit_apodization)||(h.transmit_apodization.window==uff.window.none)
                tx_apodization=ones(h.input.N_pixels,h.input.N_waves);
            else
                h.transmit_apodization.focus = h.input.scan;
                tx_apodization=h.transmit_apodization.data;
            end
            
            % building a apodization matrix
            tx_apodization=reshape(tx_apodization,[h.input.N_pixels, 1, h.input.N_waves]);
            apodization_matrix=bsxfun(@times,tx_apodization,rx_apodization);

            % compute signal phase 
            signal_phase = angle(h.input.data);
            
            % auxiliary phase
            mask=(signal_phase<=0);
            auxiliary_phase=signal_phase-pi;
            auxiliary_phase(mask)=signal_phase(mask)+pi;
            
            % declare output structure
            h.output=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data
            h.FCA=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data
            h.FCC=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data

            switch h.dimension
                case dimension.both
                    coherent_sum=sum(sum(h.input.data,2),3);
                    
                    % collapsing 2nd and 3rd dimmension into 2nd dimension
                    signal_phase=reshape(signal_phase,[h.input.N_pixels, h.input.N_channels*h.input.N_waves 1 h.input.N_frames]);
                    auxiliary_phase=reshape(auxiliary_phase,[h.input.N_pixels, h.input.N_channels*h.input.N_waves 1 h.input.N_frames]);
                    apodization_matrix=reshape(apodization_matrix,[h.input.N_pixels, h.input.N_channels*h.input.N_waves]);

                    std_phase=tools.weigthed_std(signal_phase,apodization_matrix,2);
                    std_auxiliary=tools.weigthed_std(auxiliary_phase,apodization_matrix,2);
                    std_complex=sqrt(tools.weigthed_var(cos(signal_phase),apodization_matrix,2)+tools.weigthed_var(sin(signal_phase),apodization_matrix,2));
                case dimension.transmit
                    coherent_sum=sum(h.input.data,3);
                    std_phase=tools.weigthed_std(signal_phase,apodization_matrix,3);
                    std_auxiliary=tools.weigthed_std(auxiliary_phase,apodization_matrix,3);
                    std_complex=sqrt(tools.weigthed_var(cos(signal_phase),apodization_matrix,3)+tools.weigthed_var(sin(signal_phase),apodization_matrix,3));
                case dimension.receive
                    coherent_sum=sum(h.input.data,2);
                    std_phase=tools.weigthed_std(signal_phase,apodization_matrix,2);
                    std_auxiliary=tools.weigthed_std(auxiliary_phase,apodization_matrix,2);
                    std_complex=sqrt(tools.weigthed_var(cos(signal_phase),apodization_matrix,2)+tools.weigthed_var(sin(signal_phase),apodization_matrix,2));
                otherwise
                    error('Unknown dimension mode; check HELP dimension');
            end
            % min of std deviation of phase and auxiliary phase
            sf=bsxfun(@min,std_phase,std_auxiliary);
                
            % Absolute phase coherence factor
            h.FCA.data=1-(h.gamma/h.sigma_0).*sf;
            h.FCA.data(h.FCA.data < 0) = 0;
            h.FCA.data(isnan(h.FCA.data)) = 0;

            % Complex phase coherence factor
            h.FCC.data=1-std_complex;

            % phase coherent factor image            
            h.output.data = h.FCC.data .* coherent_sum; 
            
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



