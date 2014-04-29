function varargout = fftabs2plot(data, fs, justaxis)
    % [fftData, fax] = fftabs2plot(data, fs, [false]) 
    % [fax] = fftabs2plot(data, fs, true) 
    
    % with inspiration from: 
    % http://www.mathworks.com/support/tech-notes/1700/1702.html
    % just we devide by number of samples for the scaling and no power plot
    if ~exist('justaxis','var')
        justaxis = false;
    end;
    
    nfft = length(data);
    
    % length of interesting frequency information
    halfLen = ceil((nfft+1)/2);
    
    % frequency axis
    fax = (0:halfLen-1)*fs/nfft;
        
    if justaxis
        varargout{1} = fax;
        return
    end;
        
    % do fourier transform 
    fftData = fft(data);   
    
    % keep half of it
    fftData = fftData(1:halfLen);
    
    % we're only interested in the magnitude is 
    fftData = abs(fftData);
    
    % normalize data
    fftData = fftData/nfft;
    
    varargout{1} = fftData;
    varargout{2} = fax;
end