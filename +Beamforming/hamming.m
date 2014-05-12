function w = hamming(N,n)
% HAMMING  Modified Hamming window.
%   HAmming(N) is identical to the hamming window in the signal processing toolbox,
%   and returns the N point hamming window.
%


if nargin < 2,
    n = (-(N-1)/2:(N-1)/2)';
elseif numel(n) == 1,
    n = (-(N-1)/2:(N-1)/2)' + n;
end

arg = min(1, max(0, n/(N-1)+0.5));
w = 0.54 - 0.46 * cos(2*pi*arg);
