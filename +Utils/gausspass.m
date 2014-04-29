function [y,spectrum] = gausspass(sig,fc,rel_bandwidth,fs,nonnormal)

% [y,spectrum] = gausspass(sig,fc,rel_bandwidth,fs)
% sig -- signal to be filtered. If matrix, then filter is applied to each column.
% fs -- sampling frequency
% fc -- centre frequency
% NB! rel_bandwidth -- bandwidth relative to fc!
% If fc = 0, rel_bandwidth should be replaced by absolute bandwidth

% Normalize frequencies
fc = fc/fs;
fs = 1;

if ~exist('nonnormal'),
    nonnormal = 0;
end

if fc==0,
    bandwidth = rel_bandwidth;
elseif fc > 0,
    bandwidth = fc*rel_bandwidth;
end
if fs/(fc + 0.5*bandwidth) < 2,
    warning('Sampling rate is lower than Nyquist.')
end
if size(sig,1)==1,
    sig = sig.';
end

[L,N] = size(sig);
f = (0:L-1)/L*fs;

variance = bandwidth^2/(8*log(2));

spectrum = gaussian(L,fs,fc,variance)' + gaussian(L,fs,fs-fc,variance)';
win = sqrt(spectrum);
win = win/sqrt(win'*win)*sqrt(L); % Normalize energy of window
spectrum = win.^2;

Y = fft(sig)/sqrt(L).*repmat(win,1,N);
y = ifft(Y,'symmetric')*sqrt(L); %imag(y) are almost zero
%keyboard

if nonnormal==1,
    for j = 1:N,
        y(:,j) = y(:,j) * sqrt((sig(:,j)'*sig(:,j))/(y(:,j)'*y(:,j)));
    end
end