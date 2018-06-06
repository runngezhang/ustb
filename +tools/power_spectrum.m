function [fx, F] = power_spectrum(data,fs,normalised)
    % function [fx, F] = power_spectrum(data,fs)
    
    if ~exist('normalised')
        normalised=0;
    end
    
    Nfft=size(data,1);
    F = fftshift(fft(data, Nfft));
    F = abs(squeeze(mean(mean(mean(abs(F(:,:,:,:)),4),3),2)));
    if normalised 
        F = F/max(F); 
    end
    fx = linspace(-fs/2,fs/2, Nfft);
    
end

