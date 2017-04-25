classdef phase_coherence_factor < postprocess.postprocess
    %% Phase Coherence Factor postprocessing
    properties
        gamma = 1;
    end
    
    methods
        function out_dataset = go(h,in_dataset)
            out_dataset=in_dataset(1);
            wb=waitbar(0,'Phase Coherence Factor Postprocessing');
            
            if numel(in_dataset) > 1
                data_cube = zeros(in_dataset(1).scan.N_z_axis, in_dataset(1).scan.N_x_axis,numel(in_dataset));
                for i = 1:numel(in_dataset)
                    data_cube(:,:,i) = reshape(in_dataset(i).data,[in_dataset(1).scan.N_z_axis in_dataset(1).scan.N_x_axis]);
                end
                bmf.receive_apodization.window=uff.window.none
                bmf.transmit_apodization.window=uff.window.none;
                bmf = phase_coherence_factor();
                bmf.data_cube = data_cube;
            else
                data_cube = reshape(in_dataset.data,[in_dataset(1).scan.N_z_axis in_dataset(1).scan.N_x_axis,h.bmf.channel_data.N_elements]);
                
                bmf = phase_coherence_factor(h.bmf);
                bmf.data_cube = data_cube;
                bmf.receive_apodization.probe=h.bmf.channel_data.probe;
                bmf.receive_apodization.scan=h.bmf.scan;
                rx_apo = bmf.receive_apodization.data;
                bmf.transmit_apodization.probe=h.bmf.channel_data.probe;
                bmf.transmit_apodization.scan=h.bmf.scan;
                tx_apo = bmf.transmit_apodization.data;
                
                bmf.apo = tx_apo.*rx_apo; 
            end
            bmf.gamma = h.gamma;
            
            % Call the adaptive beamformer implementation
            image = bmf.phase_coherence_factor_implementation();
            
            out_dataset.data = image;
            close(wb);
        end
    end
end
