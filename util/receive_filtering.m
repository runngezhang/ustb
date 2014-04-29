function [y,b,a]=receive_filtering(x,freq_range,fs,order)
%RECEIVE_FILTERING Filtering of the received signal
%
%[y,b,a]=receive_filtering(x,freq_range,fs,order)
%
% Function for performing receive filtering of data from PELab
%
% Input: 
%   x           - signal to be filtered
%   freq_range  - frequency range in MHz [low_freq high_freq]
%   fs          - sampling frequency of data
%   order       - filter order (Default = 2)
% Output:
%   y           - filtered signal
%   b,a         - filter coefficients

% Created: 30/8-2007 by Svein-Erik Måsøy
% Modified:

if nargin<4
    order=2;
end

Ts=1/fs;
freq_range = freq_range*1e6*Ts*2;
L=length(x);
[b,a] = butter(order,freq_range);
win = tukeywin(L,0.15);
xf = x.*win(:,ones(1,size(x,2)));
y = filtfilt(b,a,xf);
