function [xn n] = addnoise(x,SNR,fs,freq_band)

import Utils.gausspass

n = randn(size(x));

if nargin > 3    
    freq_band = freq_band/fs;
    fs = 1;
    fc = sum(freq_band)/2;
    rel_bandwidth = (freq_band(2)-freq_band(1))/fc;
    n = randn(size(x));
    for kk=1:size(x,3)
        n(:,:,kk) = gausspass(n(:,:,kk),fc,rel_bandwidth,fs);
    end
end

x0 = x/rms(x);
n0 = n/rms(n);

c = 10^(-SNR/20);

xn = x0 + c*n0;

