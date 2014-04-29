function env=envelope(fc,Nc,N,Ts,type)
% ENVELOPE
% env=envelope(fc,Nc,N,type)
% Computes the envelope
% Input:
% fc - center frequency of the pulse to transmit
% Nc - number of cycles to transmit
% N - number of samples in pulse
% Ts - sampling frequency 
% type - 'cosine', 'gaussian', or 'none'
%
% Output:
% env - the envelope
%
% Created by: Thor Andreas Tangen 2007
% Last modified by: 



switch type
    case 'cosine'
        Tp = Nc/fc;
        fe = 1/Tp;
        tenv = Ts*(0:(N-1))';
        env = 0.5-0.5*cos(2*pi*fe*tenv);        
    case 'gaussian'
        sigma = N/5;
        tenv = linspace(-N/2,N/2,N);
        env = exp(-tenv.^2'/(2*sigma^2));
    case 'none'
        env = ones(N,1);
    otherwise
        env = ones(N,1);
end