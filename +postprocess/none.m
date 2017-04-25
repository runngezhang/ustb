classdef none < postprocess.postprocess
    
    methods
        function out_dataset = go(h,in_dataset)
            %NONE Doesn't change the dataset
            
            out_dataset=in_dataset;
        end
    end
end

