function [filtered_p,h,w] = band_pass(p, Fs, F, dev)
    % function filtered_p = band_pass(p, Fs, F, dev)
    
    if ~exist('dev')
        dev=[1e-3 1e-3 1e-3];     % ripple/attenuation spec
    end
    
    % filter specification
    A=[0 1 0];                % band type: 0='stop', 1='pass'
    [M,Wn,beta,typ]= kaiserord(F,A,dev,Fs);  % window parameters
    b=fir1(M,Wn,typ,kaiser(M+1,beta),'noscale'); % filter design
    
    [h,w] = freqz(b);

    % filtering
    sz=size(p);
    filt_delay=round((length(b)-1)/2);
    filtered_p=filter(b,1,[p; zeros([filt_delay sz(2:end)])],[],1);

    % correcting the delay
    filtered_p=filtered_p((filt_delay+1):end,:,:,:);

end

