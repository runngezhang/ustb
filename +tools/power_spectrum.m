function [fx, F] = power_spectrum(data,fs)
    % function [fx, F] = power_spectrum(data,fs)
        
    Nfft=size(data,1);
    F = fftshift(fft(data, Nfft));
    F = abs(squeeze(mean(mean(mean(abs(F(:,:,:,:)),4),3),2)));
    F = F/max(F);
    fx = linspace(-fs/2,fs/2, Nfft);
    
end

