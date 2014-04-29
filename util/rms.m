function r = rms(x,w)

if nargin == 1
    r = sqrt(mean(x.*conj(x)));
elseif nargin == 2
    X = x.*conj(x);
    P = mean(w);
    M = mean(w.*X);
    r = sqrt(M/P);
end