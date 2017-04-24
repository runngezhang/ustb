classdef coherence_factor < postprocess.postprocess
    
    properties
       bmf 
    end
    
    methods
        function out_dataset = go(h,in_dataset)
            %% COHERENT_COMPOUND Coherently compounds beamformed datasets
            
            out_dataset=in_dataset(1);
            wb=waitbar(0,'Postprocessing');
            data_cube = zeros(b_data_temp(1).scan.N_z_axis, b_data_temp(1).scan.N_x_axis,numel(b_data_temp));
            for i = 1:numel(in_dataset)
                data_cube(:,:,i) = reshape(in_dataset(i).data,[in_dataset(1).scan.N_z_axis in_dataset(1).scan.N_x_axis]);
            end
            
            bmf = beamformer.coherence_factor(h.bmf);
            bmf.data_cube = data_cube;
            bmf.receive_apodization.probe=h.bmf.channel_data.probe;
            bmf.receive_apodization.scan=h.bmf.scan;
            rx_apo = bmf.receive_apodization.data;
            bmf.transmit_apodization.probe=h.bmf.channel_data.probe;
            bmf.transmit_apodization.scan=h.bmf.scan;
            tx_apo = bmf.transmit_apodization.data;
            
            bmf.apo = tx_apo.*rx_apo;
            % Call the adaptive beamformer implementation
            image = bmf.coherence_factor_implementation();
            
            out_dataset.data = image;
            close(wb);
            
        end
    end 
end
