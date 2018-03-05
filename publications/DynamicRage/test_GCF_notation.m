%Create a signal
s = [1 2 3 4 5 4 3 2 1];

% Implement with the direct DFT as written in the article
% NB, only for the "low frequency" region -M0:M0
M = length(s);
idx = 1;
M0 = 5;
P = zeros(1,length(-M0:M0));
S = zeros(1,length(-M0:M0));
tic
for a = -M0:M0
    for m = 1:M
       P(idx) =  P(idx) + exp(-j*pi*(a))*s(m)*exp(-j*2*pi*(m-1)*(a)/M);
       S(idx) =  S(idx) + s(m)*exp(-j*2*pi*(m-1-M/2)*(a)/M);
    end
    idx = idx + 1
end
toc

figure(1);clf
plot(-M0:M0,abs(P));hold on;
plot(-M0:M0,abs(S),'b*');
title(['The low frequency region M0 = ',num2str(M0)]);

tic
fft(s);
toc

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