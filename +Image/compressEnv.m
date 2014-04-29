function [Compressed] = compressEnv(aEnv, img_max, value_max, gain, dyn, method, warn, varargin)
    % [Compressed] = compressEnv(aEnv, img_max, value_max, gain, dyn,
    % method, varargin)
    % compresses the abs(aEnv) within 0:value_max to 1:img_max offering
    % different methods, usually with gain and dynamic range
    if ~exist('warn','var')
        warn = false;
    end;
    switch method
        case 'imagelog2'
            % for squared envelope data
            Compressed = imagelog2(abs(aEnv),gain,dyn,255);
        case 'imagelog'
            % for envelope data
            Compressed = img_max/dyn*(gain+20*log10(max(eps,aEnv)));
            Compressed = min(img_max, max(1,round(Compressed)));
        case 'imagelog-map'
            % faster than imagelog for big datasets
            try
                factor = 1;
                immap = img_max/dyn*(gain+20*log10(max(eps,0:1/factor:value_max)));
                immap = round(min(img_max, max(0,immap)) + 1);
                Compressed = immap(ceil(min(aEnv,value_max)*factor)+1);

                % to check compression find first value
                firstidx = max(1,find(immap>1,1)-1);
                lastidx = min(find(immap<256,1,'last'),length(immap));
                if (20*log10(lastidx/firstidx)< dyn-2) && warn
                    disp(sprintf('compresEnv: image map dynamic range (%2.2f dB) is smaller than demanded (%2.2f dB)',20*log10(lastidx/firstidx),dyn));
                end
            catch
                if ~isreal(aEnv) || ~isempty(aEnv(aEnv < 0))
                    error('aEnv has to be real and positive');
                else
                    error('An error was detected. Aborting');
                end;
            end
        otherwise
            warning('Method %s unknown',method);
    end
end
