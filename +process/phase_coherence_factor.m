classdef phase_coherence_factor < process
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
%   $Last updated: 2017/05/02$

    %% constructor
    methods (Access = public)
        function h=phase_coherence_factor()
            h.name='Phase Coherence Factor MATLAB';   
            h.reference='J. Camacho and C. Fritsch, Phase coherence imaging of grained materials, in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, vol. 58, no. 5, pp. 1006-1015, May 2011.';                
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};    
            h.version='v1.0.1';
        end
    end

    %% Additional properties
    properties
        gamma=1;                                               % mixing ratio
        sigma_0=pi/sqrt(3);                                    % reference phase value
        PCF                                                    % BEAMFORMED_DATA class with the computed phase coherent factor
    end

    methods
        function [out_data h]=go(h)
            % declare & copy beamformed dataset
            out_data=uff.beamformed_data(h.beamformed_data(1));
                        
            % initialize coherent sum
            coherent=zeros(size(out_data.data));
            
            % check if we have information about apodization
            rx_apodization=ones([size(out_data.data,1),size(h.beamformed_data,1)]);
            tx_apodization=ones([size(out_data.data,1),size(h.beamformed_data,2)]);
            if ~isempty(h.transmit_apodization)&~isempty(h.receive_apodization)&~isempty(h.channel_data.probe)
                % receive
                if size(h.beamformed_data,1)>1
                    h.receive_apodization.scan=out_data.scan;
                    h.receive_apodization.probe=h.channel_data.probe;
                    rx_apodization=h.receive_apodization.data();                
                end
                
                % transmit
                if size(h.beamformed_data,2)>1
                    h.transmit_apodization.scan=out_data.scan;
                    h.transmit_apodization.probe=h.channel_data.probe;
                    tx_apodization=h.transmit_apodization.data();
                end
            else
                warning('Missing probe and apodization data; full aperture is assumed.');
            end
            
            % start process
            tools.workbar();
            online_phase=tools.weighted_incremental_variance();
            online_auxiliary_phase=tools.weighted_incremental_variance();
            [Nrx Ntx]=size(h.beamformed_data);
            N=Nrx*Ntx;
            for n_rx=1:Nrx
                for n_tx=1:Ntx
                    % progress bar
                    n=(n_rx-1)*Ntx+n_tx;
                    if mod(n,round(N/100))==2
                        tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                    end
                    
                    % accumulate coherent signal
                    coherent=coherent+h.beamformed_data(n_rx,n_tx).data;

                    % compute signal phase
                    signal_phase = angle(h.beamformed_data(n_rx,n_tx).data);

                    % auxiliary phase
                    mask=(signal_phase<=0);
                    auxiliary_phase=signal_phase-pi;
                    auxiliary_phase(mask)=signal_phase(mask)+pi;

                    % weighted incremental variance
                    apodization=tx_apodization(:,n_tx).*rx_apodization(:,n_rx);
                    online_phase.add(apodization,signal_phase);
                    online_auxiliary_phase.add(apodization,auxiliary_phase);
                end
            end
            
            % min of std deviation of phase and auxiliary phase
            sf=bsxfun(@min,online_phase.std(),online_auxiliary_phase.std());
                
            % Phase Coherence Factor
            h.PCF = uff.beamformed_data(out_data); 
            h.PCF.data=1-(h.gamma/h.sigma_0).*sf;
            h.PCF.data(h.PCF.data < 0) = 0;
            h.PCF.data(isnan(h.PCF.data)) = 0;
                        
            % phase coherent factor image            
            out_data.data = h.PCF.data .* coherent;            
        end   
    end
end



