classdef phase_coherence_factor < postprocess.postprocess
    
    properties
       bmf 
       gamma = 1;
    end
    
    methods
        function out_dataset = go(h,in_dataset)
            %% COHERENT_COMPOUND Coherently compounds beamformed datasets
            
            out_dataset=in_dataset(1);
            wb=waitbar(0,'Postprocessing');
            for i = 1:numel(in_dataset)
                data_cube(:,:,i) = reshape(in_dataset(i).data,[in_dataset(1).scan.N_z_axis in_dataset(1).scan.N_x_axis]);
            end
            
            bmf = beamformer.phase_coherence_factor(h.bmf);
            bmf.data_cube = data_cube;
            bmf.receive_apodization.probe=h.bmf.channel_data.probe;
            bmf.receive_apodization.scan=h.bmf.scan;
            rx_apo = bmf.receive_apodization.data;
            bmf.transmit_apodization.probe=h.bmf.channel_data.probe;
            bmf.transmit_apodization.scan=h.bmf.scan;
            tx_apo = bmf.transmit_apodization.data;
            
            bmf.apo = tx_apo.*rx_apo;
            
            bmf.gamma = h.bmf.gamma;
            
            % Call the adaptive beamformer implementation
            image = bmf.phase_coherence_factor();
            
            out_dataset.data = image;
            close(wb);
            
        end
    end 
end
