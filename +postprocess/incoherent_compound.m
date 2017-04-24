classdef incoherent_compound < postprocess.postprocess
    
    methods
        function out_dataset = go(h,in_dataset)
            
            %INCOHERENT_COMPOUND Incoherently compounds beamformed datasets
            
            out_dataset=in_dataset(1);
            out_dataset.data=abs(out_dataset.data);
            wb=waitbar(0,'Postprocessing');
            for n=2:length(in_dataset)
                out_dataset.data=out_dataset.data+abs(in_dataset(n).data);
            end
            close(wb);
        end
    end
end

