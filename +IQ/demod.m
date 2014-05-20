function [iq,a,b,c,s] = demod(rf,fs,fdemod,method,varargin)
% [iq,a,b,c,s] = demod_ustb(rf,fs,fdemod,method,varargin);
% 
% IQ demodulation. There is two different ways to do this, either by hilbert transform or by using a narrowband approximation.
% For using the Hilbert transform approach, use method='hilbert'. RF is the signal to demodulate, FS is the sampling frequency 
% and FDEMOD is the demodulation frequency.
%
% IQ = DEMODULATE(RF,FS,FDEMOD,METHOD,B,A)
% The narrowband method(method='narrowband') down mixes the signal using the demodulation frequency and then filters the
% signal using the filter specified by B and A.
%
% [IQ,A,B] = DEMODULATE(RF,FS,FDEMOD,METHOD,NTAPS,CUTOFF)
% The narrowband method(method='narrowband') down mixes the signal using the demodulation frequency, computes a FIR low pass
% filter using the specified order and cutoff frequency in MHz and the filters the signal with this filter. The filter is computed
% using fir1 and a hamming window.
reshaped = false;

if ndims(rf)>2
    saveshape = size(rf);
    rf = reshape(rf,size(rf,1),[]);
    reshaped = true;
end;

switch method
    case 'hilbert'
        if nargin == 4
            fs_iq = fs;
        end
        
        I = 1;%max([round(fs/fs_iq) 1]);

        [N,M]=size(rf);        

        if length(fdemod) > 1
            fdemod = linspace(fdemod(1),fdemod(2),N)';
        end
        
        t = (0:(size(rf,1)-1))'/fs;
        mix = exp(-1i*2*pi*fdemod.*t);        
        mix = mix(:,ones(1,M));
        
        dm = hilbert(rf).*mix;
        iq = dm(1:I:end,:);
        
        a = 1;
        b = 1;
        a = mix;
    case 'narrowband'
        if nargin < 5
            error('');
        end
        
        if ndims(rf) > 1
            [L M] = size(rf);
        else
            L = length(rf);
            M = 1;
            rf = reshape(rf,L,1);
        end
        
        if length(fdemod) > 1
            fdemod = linspace(fdemod(1),fdemod(2),L)';
        end
        
        Ts = 1/fs;
        t = Ts*(0:(L-1))';
        c = cos(-2*pi*fdemod.*t);        
        s = sin(-2*pi*fdemod.*t);                
               
%         NN = 50;
%         phi = linspace(-pi*fdemod*Ts,0,NN);
%         for kk=1:NN            
%             c_kk = cos(-2*pi*fdemod*t + phi(kk));                        
%             c2_kk = sign(c_kk).*2.^round(log2(abs(c_kk)));        
%             err(kk) = mean((c_kk - c2_kk).^2);
%         end
%         [dummy ix] = min(err);
        
%         c = cos(-2*pi*fdemod*t + 0*phi(ix));
%         s = sin(-2*pi*fdemod*t + 0*phi(ix));
%         c = sign(c).*2.^round(log2(abs(c)));        
%         s = sign(s).*2.^round(log2(abs(s)));        
        
%         figure(100)
%         plot(1:L,cc,1:L,sign(cc).*2.^round(log2(abs(cc))));
%         figure(101)
%         plot(phi/pi/fdemod/Ts/2,err)
%         figure(102)
%         plot(round(log2(abs(cc))))
%         
        c = c(:,ones(1,M));        
        s = s(:,ones(1,M));     
        
        re = rf.*c;
        im = rf.*s;
        
        if nargin >= 6
            if length(varargin{1}) > 1
                b = varargin{1};
                a = varargin{2};
            else
                Ntaps = varargin{1};
                cutoff = varargin{2};
                cutoff = cutoff*Ts*2;
                b = fir1(Ntaps,cutoff,hamming(Ntaps+1),'scale');
                a = 1;                              
            end
        end
        
        %b = b/sqrt(b*b');
        %a = a/sqrt(a*a');        
        
        re = 0.5*filter(b,a,re(:));            
        im = 0.5*filter(b,a,im(:));
        re = reshape(re,size(rf,1),size(rf,2));
        im = reshape(im,size(rf,1),size(rf,2));
        
        iq = re + 1i*im;
end
if reshaped
    iq = reshape(iq, saveshape);
end;
