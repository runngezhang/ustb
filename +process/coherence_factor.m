classdef coherence_factor < process
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
        function h=coherence_factor()
            h.name='Coherence Factor MATLAB';   
            h.reference= 'R. Mallart and M. Fink, Adaptive focusing in scattering media through sound-speed inhomogeneities: The van Cittert Zernike approach and focusing criterion, J. Acoust. Soc. Am., vol. 96, no. 6, pp. 3721-3732, 1994';                
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};    
            h.version='v1.0.2';
        end
    end
        
    %% Additional properties
    properties
        CF                                            % BEAMFORMED_DATA class with the computed coherent factor
        active_element_criterium=0.16;                % value to decide whether an element is used or not. This value depends on the SNR so it must be adjusted on a case-by-case basis.
        operation = operation.both;                   % OPERATION class that specifies whether the process will run only on transmit, receive, or both.
    end

    methods
        function [out_data h]=go(h)
           
            % check if we have information about apodization
            rx_apodization=ones([h.beamformed_data(1).N_pixels,size(h.beamformed_data,1)]);
            tx_apodization=ones([h.beamformed_data(1).N_pixels,size(h.beamformed_data,2)]);
            if ~isempty(h.transmit_apodization)&~isempty(h.receive_apodization)&~isempty(h.channel_data.probe)
                % receive
                if size(h.beamformed_data,1)>1
                    h.receive_apodization.scan=h.beamformed_data(1).scan;
                    h.receive_apodization.probe=h.channel_data.probe;
                    rx_apodization=h.receive_apodization.data();                
                end
                
                % transmit
                if size(h.beamformed_data,2)>1
                    h.transmit_apodization.scan=h.beamformed_data(1).scan;
                    h.transmit_apodization.probe=h.channel_data.probe;
                    tx_apodization=h.transmit_apodization.data();
                end
            else
                warning('Missing probe and apodization data; full aperture is assumed.');
            end
            
            % declare temporary variables            
            [Nrx Ntx]=size(h.beamformed_data); 
            switch h.operation
                case operation.both
                    coherent=zeros(h.beamformed_data(1).N_pixels,1,1);
                    incoherent=zeros(h.beamformed_data(1).N_pixels,1,1);
                    M=zeros(h.beamformed_data(1).N_pixels,1,1);
                case operation.transmit
                    coherent=zeros(h.beamformed_data(1).N_pixels,Nrx,1);
                    incoherent=zeros(h.beamformed_data(1).N_pixels,Nrx,1);
                    M=zeros(h.beamformed_data(1).N_pixels,Nrx,1);
                case operation.receive
                    coherent=zeros(h.beamformed_data(1).N_pixels,1,Ntx);
                    incoherent=zeros(h.beamformed_data(1).N_pixels,1,Ntx);
                    M=zeros(h.beamformed_data(1).N_pixels,1,Ntx);
                otherwise
                    error('Unknown operation mode; check HELP OPERATION');
            end
            
            % loop
            n=1; N=Ntx*Nrx; tools.workbar(); 
            for n_rx=1:Nrx
                for n_tx=1:Ntx
                    % progress bar
                    if mod(n,round(N/100))==1
                        tools.workbar(n/N,sprintf('%s (%s)',h.name,h.version),'USTB');
                    end
                    n=n+1;
        
                    switch h.operation
                        case operation.both
                            coherent=coherent+h.beamformed_data(n_rx,n_tx).data;
                            incoherent=incoherent+abs(h.beamformed_data(n_rx,n_tx).data).^2;
                            M=M+double((tx_apodization(:,n_tx).*rx_apodization(:,n_rx))>h.active_element_criterium);
                        case operation.transmit
                            coherent(:,n_rx)=coherent(:,n_rx)+h.beamformed_data(n_rx,n_tx).data;
                            incoherent(:,n_rx)=incoherent(:,n_rx)+abs(h.beamformed_data(n_rx,n_tx).data).^2;
                            M(:,n_rx)=M(:,n_rx)+double((tx_apodization(:,n_tx).*rx_apodization(:,n_rx))>h.active_element_criterium);
                        case operation.receive
                            coherent(:,1,n_tx)=coherent(:,1,n_tx)+h.beamformed_data(n_rx,n_tx).data;
                            incoherent(:,1,n_tx)=incoherent(:,1,n_tx)+abs(h.beamformed_data(n_rx,n_tx).data).^2;
                            M(:,1,n_tx)=M(:,1,n_tx)+double((tx_apodization(:,n_tx).*rx_apodization(:,n_rx))>h.active_element_criterium);
                    end
                            
                end
            end
            tools.workbar(1);
            
            switch h.operation
                case operation.both
                    % declare variables
                    out_data = uff.beamformed_data(h.beamformed_data(1)); 
                    h.CF = uff.beamformed_data(h.beamformed_data(1));
                    
                    % Coherent Factor
                    h.CF.data = abs(coherent).^2./incoherent./M; 
                    % coherent factor image            
                    out_data.data = h.CF.data .* coherent;
                    
                case operation.receive
                    % transmit loop
                    for ntx=1:Ntx
                        out_data(1,ntx) = uff.beamformed_data(h.beamformed_data(1)); 
                        CF(1,ntx) = uff.beamformed_data(h.beamformed_data(1));
                    
                        % Coherent Factor
                        CF(1,ntx).data = abs(coherent(:,1,ntx)).^2./incoherent(:,1,ntx)./M(:,1,ntx); 
                        
                        % coherent factor image            
                        out_data(1,ntx).data = CF(1,ntx).data .* coherent(:,1,ntx);
                    end
                    h.CF=CF;
                case operation.transmit
                    % receive loop
                    for nrx=1:Nrx
                        out_data(nrx) = uff.beamformed_data(h.beamformed_data(1)); 
                        CF(nrx) = uff.beamformed_data(h.beamformed_data(1));
                    
                        % Coherent Factor
                        CF(nrx).data = abs(coherent(:,nrx)).^2./incoherent(:,nrx)./M(:,nrx); 
                        
                        % coherent factor image            
                        out_data(nrx).data = CF(nrx).data .* coherent(:,nrx);
                    end                    
                    h.CF=CF;
                otherwise
                    error('Unknown operation type; check HELP OPERATION');
            end
        end   
    end
end



