function demodulate(h,plot_on,modulation_frequency,bandpass_frequency_vector,downsample_frequency)
  %DEMODULATE    Demodulates RF data in a dataset class
%
%   Syntax:
%   demodulate(plot_on,modulation_frequency,bandpass_frequency_vector)
%       plot_on                     Flag to plot signal power spectrum at all steps (Optional)
%       modulation_frequency        Modulation frequency in Hz. (Optional)
%       bandpass_frequency_vector   Vector with the bandpass frequencies in Hz [low_f_off low_f_on high_f_on high_f_off] (optional)
%       downsample_frequency        Sample frequency of the demodulated signal in Hz. (Optional) 
%
%   See also DATASET
    
    assert(h.format==E.signal_format.RF,'Signal format is not RF. Only RF data can be demodulated.');

    if ~exist('plot_on') plot_on=false; end
    
    %% computing central frequency and bandwidth
    disp('Estimating power spectrum');
    [fx pw]=tools.power_spectrum(h.data,h.sampling_frequency);
    fpw=filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)];
    [dc ic]=max(fpw.*(fx>0).'); fc=fx(ic);
    bw_up=min(fx((fx>fc)&(fpw<dc/2).')); % -6dB upper limit
    bw_do=max(fx((fx<fc)&(fpw<dc/2).')); % -6dB down limit
    fc=(bw_up+bw_do)/2;                  % center frequency
    bw=2*(bw_up-fc);

    if ~exist('modulation_frequency') modulation_frequency=fc; end
    
    if(plot_on)
        figure();
        subplot(1,2,1);
        plot(fx,pw,'k'); hold on; axis manual;
        axis([-2*bw_up 2*bw_up 0 1]);
    end

    %% band pass filtering
    disp('Band Pass filtering');
    if ~exist('bandpass_frequency_vector') 
        transition=bw/10;
        low_freq=max([0 fc-2*bw]);      
        high_freq=min([h.sampling_frequency/2*0.99 fc+2*bw]);
        bandpass_frequency_vector=[low_freq low_freq+transition high_freq-transition high_freq];
    end
    data = tools.band_pass(h.data,h.sampling_frequency,bandpass_frequency_vector);

    if(plot_on)
        [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
        plot([-fliplr(bandpass_frequency_vector) bandpass_frequency_vector],[0 1 1 0 0 1 1 0],'b');
        plot(fx,pw,'r--');
        plot([modulation_frequency modulation_frequency],[0 1],'m:');
        legend('Original signal','Bandpass filter','Filtered signal','Modulation frequency');
        title('Bandpass processing')
    end

    %% demodulation
    mod_sig=exp(-j*2*pi*modulation_frequency*h.time)*ones(1,size(data,2)); % demodulation sognal
    wb = waitbar(0, 'Demodulating...');
    for f=1:size(data,4)
        for n=1:size(data,3)
            waitbar(n / size(data,3), wb);
            data(:,:,n,f)=data(:,:,n,f).*mod_sig;
        end
    end
    close(wb);

    if(plot_on)    
        [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
        subplot(1,2,2);
        plot(fx,pw,'k-'); hold on; axis manual;
        axis([-2*bw_up 2*bw_up 0 1]);
    end

    %% low pass filtering
    disp('Base Band filtering');
    baseband_frequency_vector=[0.9*modulation_frequency modulation_frequency];
    data=tools.low_pass(data,h.sampling_frequency,baseband_frequency_vector);

    if(plot_on)    
        [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
        plot([-fliplr(baseband_frequency_vector) baseband_frequency_vector],[0 1 1 0],'b');
        plot(fx,pw,'r--'); 
    end

    %% resampling
    if ~exist('downsample_frequency')
        downsample_frequency=round(4*bw/1e6)*1e6; % new sampling frequency
    end
    
    dt=1/downsample_frequency;
    t=(h.time(1):dt:h.time(end));
    data=reshape(data,size(h.data)); % we get back singleton dimmensions
    downsample_data=zeros(length(t),size(data,2),size(data,3),size(data,4));
    wb = waitbar(0, 'Resampling...');
    for f=1:size(data,4)
        for ntx=1:size(data,3)
            for nrx=1:size(data,2)
                downsample_data(:,nrx,ntx,f)=interp1(h.time,data(:,nrx,ntx,f),t,'linear',0);
            end
            waitbar(ntx / size(data,3), wb);
        end
    end
    close(wb);

    if(plot_on)  
        [fx pw] = tools.power_spectrum(downsample_data,downsample_frequency);
        plot(fx,pw,'g:'); 
        legend('Demodulated signal','Baseband filter','Filtered signal','Resampled signal');
        title('Baseband processing')
    end

    %% write in the output 
    h.format=E.signal_format.IQ;
    h.modulation_frequency=modulation_frequency;
    h.time=t.';
    h.data=downsample_data;
end








