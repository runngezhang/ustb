classdef phase_coherence_factor_alternative < process
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
    %   both. Use dimension do decide wihich.
    %
    %   implementers: Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %                 Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
    %
    %   $Last updated: 2017/05/05$
    
    %% constructor
    methods (Access = public)
        function h=phase_coherence_factor_alternative()
            h.name='Phase Coherence Factor MATLAB';
            h.reference='J. Camacho and C. Fritsch, Phase coherence imaging of grained materials, in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, vol. 58, no. 5, pp. 1006-1015, May 2011.';
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>','Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};
            h.version='v1.0.2';
        end
    end
    
    %% Additional properties
    properties
        gamma=1;                                               % mixing ratio
        sigma_0=pi/sqrt(3);                                    % reference phase value
        PCF                                                    % BEAMFORMED_DATA class with the computed phase coherent factor
        dimension = dimension.both;                            % dimension class that specifies whether the process will run only on transmit, receive, or both.
    end
    
    methods
        function [out_data h]=go(h)
            assert(~isempty(h.dimension),'Please select dimension to process, for MV transmit and receive is deinfed.');
            
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
                    h.transmit_apodization.sequence = h.channel_data.sequence;
                    h.transmit_apodization.scan=h.beamformed_data(1).scan;
                    h.transmit_apodization.probe=h.channel_data.probe;
                    tx_apodization=h.transmit_apodization.data();
                end
            else
                warning('Missing probe and apodization data; full aperture is assumed.');
            end
            
            tools.workbar();
            [Nrx Ntx]=size(h.beamformed_data);
            tools.workbar(1);
            %% Call the phase coherence implementation
            switch h.dimension
                case dimension.both
                    % declare & copy beamformed dataset
                    out_data=uff.beamformed_data(h.beamformed_data(1));    
                    % initialize coherent sum
                    coherent=zeros(size(out_data.data));
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
                    tools.workbar(1);
                    
                    % min of std deviation of phase and auxiliary phase
                    sf=bsxfun(@min,online_phase.std(),online_auxiliary_phase.std());
                    
                    % Phase Coherence Factor
                    h.PCF = uff.beamformed_data(out_data);
                    h.PCF.data=1-(h.gamma/h.sigma_0).*sf;
                    h.PCF.data(h.PCF.data < 0) = 0;
                    h.PCF.data(isnan(h.PCF.data)) = 0;
                    
                    % phase coherent factor image
                    out_data.data = h.PCF.data .* coherent;
                case dimension.transmit
                    %%
                    for n_rx = 1:Nrx
                        % Create apodization matrix to know which elements
                        % are "active"
                        apo_matrix = reshape(bsxfun(@times,tx_apodization,rx_apodization(:,n_rx)),h.beamformed_data(1).scan.N_z_axis,h.beamformed_data(1).scan.N_x_axis,Ntx);
                        % Pick out correct data_cube
                        for n_tx = 1:Ntx
                            data_cube(:,:,n_tx) = (reshape(h.beamformed_data(n_rx,n_tx).data,h.beamformed_data(1).scan.N_z_axis,h.beamformed_data(1).scan.N_x_axis));
                        end
                        
                        %A hack to set non active elements to zero for the
                        %alpinion scanner FI who only use 64 active
                        %elements
                        if isempty(h.channel_data.N_active_elements) && sum(h.channel_data.N_active_elements ~= h.channel_data.N_elements)
                            apo_matrix(abs(data_cube)<eps) = 0;
                        end
                        
                        % Calculate PCF
                        pcf = h.phase_coherence_implementation(data_cube,apo_matrix,[num2str(n_rx),'/',num2str(Nrx)]);
                        % Create coherent DAS image
                        das = sum(data_cube,3);
                        % Mulitply DAS and PCF
                        image = das.*pcf;
                        % declare & copy beamformed dataset
                        out_data(n_rx)=uff.beamformed_data(h.beamformed_data(1));
                        out_data(n_rx).data = image(:);
                        PCF(1,n_tx) = uff.beamformed_data(h.beamformed_data(1));
                        PCF(1,n_tx).data = pcf(:);
                    end
                    h.PCF = PCF;
                case dimension.receive
                    for n_tx = 1:Ntx
                        % Create apodization matrix to know which elements
                        % are "active"
                        apo_matrix = reshape(bsxfun(@times,tx_apodization(:,n_tx),rx_apodization),h.beamformed_data(1).scan.N_z_axis,h.beamformed_data(1).scan.N_x_axis,Nrx);
                        
                        % Pick out correct data_cube
                        for n_rx = 1:Nrx
                            data_cube(:,:,n_rx) = (reshape(h.beamformed_data(n_rx,n_tx).data,h.beamformed_data(1).scan.N_z_axis,h.beamformed_data(1).scan.N_x_axis));
                        end
                        
                        %A hack to set non active elements to zero for the
                        %alpinion scanner FI who only use 64 active
                        %elements
                        if isempty(h.channel_data.N_active_elements) && sum(h.channel_data.N_active_elements ~= h.channel_data.N_elements)
                            apo_matrix(abs(data_cube)<eps) = 0;
                        end
                        
                        % Calculate PCF
                        pcf = h.phase_coherence_implementation(data_cube,apo_matrix,[num2str(n_tx),'/',num2str(Ntx)]);
                        % Create coherent DAS image
                        das = sum(data_cube,3);
                        
                        % Mulitply DAS and PCF
                        image = das.*pcf;
                        % declare & copy beamformed dataset
                        out_data(1,n_tx)=uff.beamformed_data(h.beamformed_data(1));
                        out_data(1,n_tx).data = image(:);
                        PCF(1,n_tx) = uff.beamformed_data(h.beamformed_data(1));
                        PCF(1,n_tx).data = pcf(:);
                    end
                    h.PCF = PCF;
                otherwise
                    error('Unknown dimension mode; check HELP dimension');
            end
            
            
        end
        
        % The actual implementation of the phase coherence working on a
        % "data cube"
        function pcf = phase_coherence_implementation(h,data_cube,apod_matrix,progress)
            Zs = size(data_cube,1);
            tools.workbar();
            sf = zeros(size(data_cube,1),size(data_cube,2));
            for zs = 1:Zs
                tools.workbar(zs/Zs,sprintf('%s %s (%s)',h.name,h.version,progress),'Phase Coherence Factor');
                for xs = 1:size(data_cube,2)
                    phase=angle(data_cube(zs,xs,logical(squeeze(apod_matrix(zs,xs,:)))));
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
            tools.workbar(1);
            % APCF
            pcf=1-sf/h.sigma_0;
        end
    end
end



