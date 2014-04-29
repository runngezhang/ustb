function [rf I] = iq2rf(iq,fs_iq,fdemod,fs_rf)
% [iq I] = IQ2RF(rf,fs_rf,fdemod,fs_iq);
% RF signal reconstruction
% Input:
% rf - rf signal
% fs_iq - sampling frequency of iq signal
% fdemod - demodulation frequency
% fs_rf - (optional) desired sampling frequency of rf signal. If not
% specified no upsampling is performed. If higher than fs_rf it is set equal to
% fs_iq

if nargin == 3
    fs_rf = fs_iq;
end

I = max([round(fs_rf/fs_iq) 1]);

order=3;
cutoff=1;

[N,M] = size(iq);
rf = zeros(N*I,M);

t = (0:N-1)'/fs_rf;
mix = exp(i*2*pi*fdemod*t);

for n=1:M,
   rf(:,n)=sqrt(2)*real(interp(iq(:,n),I,order,cutoff).*mix);
end