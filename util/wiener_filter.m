function G=wiener_filter(X,Y,Ns,method)
% WIENER_FILTER
% Computes a wiener filter based on the 2 spectrums.
% G = wiener_filter(X,Y,Ns,method)
% G = (X/Y)*(1/(1 + Ns/|Y|^2))
% Input:
% X - Denominator spectrum
% Y - Numerator spectrum
% Ns -  Noise parameter
% method - (optional) 'abs' or 'cplx'

if nargin == 3
    method = 'abs';
end

K = 1./(1 + Ns./(Y.*conj(Y)));
switch method
    case 'abs'
        %G = K.*(abs(X)./abs(Y));
        G = K.*(abs(X./(Y+eps)));
    case 'cplx'
        G = K.*(X./(Y+eps));
    otherwise
        G = K.*(abs(X)./(abs(Y)+eps));
end