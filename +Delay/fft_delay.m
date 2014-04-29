function varargout=fft_delay(varargin)
% FFT_DELAY
% Estimates delay based on frequency spectrums. Uses a weighted least
% squares fit to the phase difference of the spectrums. The weight is the
% multiplicated absolute values of the spectrums.
%
% tau=fft_delay(X1,X2,f,flim,unit)
% Input
% X1 - frequency spectrum of first signal
% X2 - frequency spectrum of second signal
% f - frequency axis
% lim - (optional) Specify interval to be used for the least
% squares fit.
% unit - (optional) unit of the limit, 'Hz' of 'dB'
%
% y = fft_delay(x,tau,fs)
% Input:
% x - signal to be delayed
% tau - delay
% fs - optional, if not specfied, unit sampling is assumed and tau is in
% assumed to be in fractions of samples
if nargin < 3 || (nargin == 3 && length(varargin{3}) == 1)
    % Do fequency delaying:
    x = varargin{1};
    tau = varargin{2};
    
    if nargin == 3
        fs = varargin{3};
    else
        fs = 1;
    end
    
    N = length(x);
    
    f = fs*[0:floor(N/2),-ceil(N/2)+1:-1]/N;
    
    y = real(ifft(fft(x).*exp(i*2*pi*f*tau)));
    
    if nargout == 1
        varargout{1} = y;
    end
else
    if nargin < 3
        error('Not enough inputs');
    end
    
    X1 = varargin{1};
    X2 = varargin{2};
    f = varargin{3};
    
    if nargin < 4
        lim = [f(1) f(end)];
    else
        lim = varargin{4};
    end
    
    if nargin < 5
        unit = 'Hz';
    else
        unit = varargin{5};
    end

    
    [N M] = size(X1);
    switch unit
        case 'Hz'
            ix = find_indx(f,lim);

            X1 = X1(ix(1):ix(2),:);
            X2 = X2(ix(1):ix(2),:);
            f = f(ix(1):ix(2));
            w = 2*pi*f;

            for m=1:M
                weight = abs(X1(:,m)).*abs(X2(:,m));
                weight = weight/sum(weight);
                
                phi = angle(X1(:,m).*conj(X2(:,m)));                
                tau(m) = -sum(phi.*weight.*w)./sum(weight.*w.*w);
            end
        case 'dB'
            for m=1:M
                ix = find_dB_limits(X1(1:floor(size(X1,1)/2),m),lim);

                X1m = X1(ix(1):ix(2),m);
                X2m = X2(ix(1):ix(2),m);
                w = 2*pi*f(ix(1):ix(2));

                weight = abs(X1(ix(1):ix(2),m)).*abs(X2(ix(1):ix(2),m));
                weight = weight/sum(weight);
                
                phi = angle(X1m.*conj(X2m));

                tau(m) = -sum(phi.*weight.*w)/sum(weight.*w.*w);
            end
    end
    
    if nargout == 1
        varargout{1} = tau;
    end
end