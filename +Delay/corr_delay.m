function [tau m_ix R] = corr_delay(rf1,rf2,Ts,maxlag,avg)
% CORR_DELAY
% tau = corr_delay(rf1,rf2,Ts,mlag)
% Computes the delay using cross correlation. Uses parabolic interpolation
% around the peak of the correlation function to find subsample delay.
% Input:
% rf1 - signal 1
% rf2 - signal 2
% Ts - 1/sampling frequency
% mlag - (optional) maximum expected delay in samples. Speeds up
% computation. Same as MAXLAG in xcorr.

if nargin < 4
    maxlag = max([length(rf1) length(rf2)]);
end

if nargin < 5
    avg = 0;
end

[N M] = size(rf1);
% if N < M
%     rf1 = rf1';
%     rf2 = rf2';    
% end

if size(rf1) ~= size(rf2)
    error('Signals must have same size');
end

if avg
    RF1 = fft(rf1,[],1);
    RF2 = fft(rf2,[],1);
    
    R = ifft(mean(RF1.*conj(RF2),2),[],1);
    
    if isreal(rf1)
        R = real(R);    
    end
    
    R = [R(end-maxlag+1:end,:);R(1:maxlag+1,:)];
    lags = -maxlag:maxlag;
    
    [dummy m_ix] = max(R);
    
    if m_ix <= 1 || m_ix >= (length(R)-1)
        tau = 0;
    else
        tau = Ts*((R(m_ix-1) - R(m_ix+1))/(2*(R(m_ix-1)-2*R(m_ix)+R(m_ix+1))) + lags(m_ix));   
    end        
else
    for m=1:M
        [R lags] = xcorr(rf1(:,m),rf2(:,m),maxlag);    
        [dummy m_ix] = max(R);
        if m_ix <= 1 || m_ix >= (length(R)-1)
            tau(m) = 0;
        else
            tau(m) = Ts*((R(m_ix-1) - R(m_ix+1))/(2*(R(m_ix-1)-2*R(m_ix)+R(m_ix+1))) + lags(m_ix));   
        end
    end
end
    

