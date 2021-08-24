classdef image_postprocessing < postprocess
    %   Adaptive Receive Side Compounding
    %   
    %
    %   Authors: Are Charles Jense <arej@ifi.uio.no>, Andreas Austeng <AndreasAusteng@ifi.uio.no>
    %
    %   $Last updated: 2020/17/01$
    
    %% constructor
    methods (Access = public)
        function h=image_postprocessing()
            h.name='Image Postprocessing';
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>'};
            h.version='v1.0.0';
        end
    end
    
    %% Additional properties
    properties
        dimension = dimension.receive;                % dimension class that specifies whether the process will run only on transmit, receive, or both.
        channel_data                                  % Channel data 
        scan
        median_m = 5;
        median_n = 5;
    end
    
    methods
        function [output]=go(h)       
            % check if we can skip calculation
            if h.check_hash()
                output = h.output; 
                return;
            end  
            
            % check dimensions
            assert(size(h.input.data,2)==1,'image_postprocessing only works on combined images');
            assert(size(h.input.data,3)==1,'image_postprocessing only works on combined images');

            % declare output structure
            h.output=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data
            
            % temp datastructure
            aux_data=zeros(h.input.N_pixels,1,h.input.N_waves,h.input.N_frames);
            for n_frame = 1:h.input.N_frames
                % Reshape to image matrix
                if isa(h.scan,'uff.linear_scan')
                    img = reshape(h.input.data(:,:,:,n_frame),h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_channels);
                    
                elseif isa(h.scan,'uff.sector_scan')
                    img = reshape(h.input.data(:,:,:,n_frame),h.input(1).scan.N_depth_axis,h.input(1).scan.N_azimuth_axis,h.input.N_channels);
                    
                end
                image = image_postprocessing_implementation(h,img,[num2str(n_frame),'/',num2str(h.input.N_frames)]);
                aux_data(:,1,1,n_frame) = image(:);
                
            end
            h.output.data = aux_data;
            
            % pass reference
            output = h.output;
            
            % update hash
            h.save_hash();
        end
        
        function img_out = image_postprocessing_implementation(h,img_in,progress)
            [N,E,M] = size(img_in);
            
            envelope_img = abs(img_in);
            
            % As an example, simply just do median filtering on the image
            img_out = medfilt2(envelope_img,[h.median_m h.median_n]);
        end

        
    end
end



