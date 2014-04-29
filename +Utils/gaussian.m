function w = gaussian(N,L,mn,vrnc)

%Creates N samples of gaussian window of length L, mean mn and variance vrnc.

k = (0:N-1)/N * L;
w = 1/(sqrt(2*pi*vrnc)) * exp(-(k-mn).^2/(2*vrnc));