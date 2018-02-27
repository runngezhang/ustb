%Create a signal
s = [1 2 3 4 5 4 3 2 1];

% Implement with the direct DFT as written in the article
% NB, only for the "low frequency" region -M0:M0
M = length(s);
idx = 1;
M0 = 2;
p = zeros(1,length(-M0:M0));
for a = -M0:M0
    for m = 1:M
       p(idx) =  p(idx) + exp(-j*pi*(a))*s(m)*exp(-j*2*pi*(m-1)*(a)/M);
    end
    idx = idx + 1
end

figure(1);clf
plot(-M0:M0,abs(p));
title(['The low frequency region M0 = ',num2str(M0)]);


%Plot the FFT and the shifted FFT
figure(2);clf
subplot(211);hold all;
plot(abs(fft(s)));
plot([M0+1 M0+1],[0 20],'r');
plot([length(s)+1-M0 length(s)+1-M0],[0 20],'r');
title('FFT');
subplot(212)
plot(abs(fftshift(fft(s))))
title('fftshifted FFT');