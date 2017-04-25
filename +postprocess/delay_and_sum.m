classdef delay_and_sum < postprocess.postprocess
    
    properties
    end
    
    methods
        function out_dataset = go(h,in_dataset)
            %% Sum a delayed inter_dataset 
            
            out_dataset = in_dataset(1);
            out_dataset.data = sum(in_dataset.data,2);
            
        end
    end 
end
