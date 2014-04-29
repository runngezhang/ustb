function xn = unit_gain(x,fc,fs)

% Make filter unit gain at center frequency
f = get_freqs(length(x),1/fs);
normval = interp1(f',abs(fft(x)),fc,'spline',0);
xn = x/normval;