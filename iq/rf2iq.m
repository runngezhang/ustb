function [iq I] = rf2iq(rf,fs_rf,fdemod,fs_iq)
% [iq I] = RF2IQ(rf,fs_rf,fdemod,fs_iq);
% IQ demodulation.
% Input:
% rf - rf signal
% fs_rf - sampling frequency of rf signal
% fdemod - demodulation frequency
% fs_iq - (optional) desired sampling frequency of iq signal. If not
% specified not decimation is used. If higher than fs_rf it is set equal to
% fs_rf

if nargin == 3
    fs_iq = fs_rf;
end

I = max([round(fs_rf/fs_iq) 1]);

[N,M]=size(rf);
iq = zeros(N/I,M);

t = (0:(size(rf,1)-1))'/fs_rf;
mix = exp(-i*2*pi*fdemod*t);

for n=1:M,
   dm = hilbert(rf(:,n)).*mix;
   iq(:,n) = dm(1:I:end);
end;