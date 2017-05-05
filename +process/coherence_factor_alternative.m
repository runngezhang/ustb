classdef coherence_factor_alternative < process
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
    %   $Last updated: 2017/05/02$
    
    %% constructor
    methods (Access = public)
        function h=coherence_factor_alternative()
            h.name='Coherence Factor MATLAB';
            h.reference= 'R. Mallart and M. Fink, Adaptive focusing in scattering media through sound-speed inhomogeneities: The van Cittert Zernike approach and focusing criterion, J. Acoust. Soc. Am., vol. 96, no. 6, pp. 3721-3732, 1994';
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};
            h.version='v1.0.2';
        end
    end
    
    %% Additional properties
    properties
        CF                                            % BEAMFORMED_DATA class with the computed coherent factor
        active_element_criterium=0.16;                % value to decide whether an element is used or not
    end
    
    methods
        function [out_data h]=go(h)
            % declare & copy beamformed dataset
            out_data=uff.beamformed_data(h.beamformed_data(1));
            
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
            
            tools.workbar();
            M=zeros(size(out_data.data));
            [Nrx Ntx]=size(h.beamformed_data);
            
            if Nrx > 1 %Do CF "on receive"
                N=Nrx*Ntx;
                for n_tx=1:Ntx
                    % initialize coherent & coherent sum
                    coherent=zeros(size(h.beamformed_data(1).data));
                    incoherent=zeros(size(abs(h.beamformed_data(1).data)));
                    M = 0;
                    for n_rx=1:Nrx
                        % progress bar
                        n=(n_tx-1)*Nrx+n_rx;
                        if mod(n,round(N/100))==2
                            tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                        end
                        
                        coherent=coherent+h.beamformed_data(n_rx,n_tx).data;
                        incoherent=incoherent+abs(h.beamformed_data(n_rx,n_tx).data).^2;
                        M=M+double((tx_apodization(:,n_tx).*rx_apodization(:,n_rx))>h.active_element_criterium);
                        
                    end
                    % declare & copy beamformed dataset
                    out_data(1,n_tx)=uff.beamformed_data(h.beamformed_data(1));  
                    
                    % Coherent Factor
                    h.CF{1,n_tx} = uff.beamformed_data();
                    h.CF{1,n_tx}.scan = out_data(1,n_tx).scan;
                    h.CF{1,n_tx}.data = abs(coherent).^2./incoherent./M;
                    
                    % coherent factor image
                    out_data(1,n_tx).data = h.CF{1,n_tx}.data .* coherent;
                end
                tools.workbar(1);
            elseif Nrx == 1 && Ntx > 1 %Do CF "on transmit" / on "compounding"
                N=Nrx*Ntx;
                % initialize coherent & coherent sum
                coherent=zeros(size(h.beamformed_data(1).data));
                incoherent=zeros(size(abs(h.beamformed_data(1).data)));
                M = 0;
                for n_tx=1:Ntx
                    for n_rx=1:Nrx
                        % progress bar
                        n=(n_rx-1)*Ntx+n_tx;
                        if mod(n,round(N/100))==2
                            tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                        end
                        
                        coherent=coherent+h.beamformed_data(n_rx,n_tx).data;
                        incoherent=incoherent+abs(h.beamformed_data(n_rx,n_tx).data).^2;
                        M=M+double((tx_apodization(:,n_tx).*rx_apodization(:,n_rx))>h.active_element_criterium);
                        
                    end
                end
                tools.workbar(1);
                % declare & copy beamformed dataset
                out_data =uff.beamformed_data(h.beamformed_data(1));
                
                % Coherent Factor
                h.CF = uff.beamformed_data(out_data);
                h.CF.data = abs(coherent).^2./incoherent./M;
                
                % coherent factor image
                out_data.data = h.CF.data .* coherent;
                
            end
        end
    end
end



