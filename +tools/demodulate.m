function out=demodulate(s)
    import tools.* 

    Fs=1/mean(diff(s.time));            % sampling frequency

    % computing central frequency and bandwidth
    disp('Estimating central frequency and bandwidth');
    [fx pw]=tools.power_spectrum(s.data,Fs);
    fpw=filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)];
    [dc ic]=max(fpw.*(fx>0).'); fc=fx(ic);
    bw_up=min(fx((fx>fc)&(fpw<dc/2).')); % -6dB upper limit
    bw_do=max(fx((fx<fc)&(fpw<dc/2).')); % -6dB down limit
    fc=(bw_up+bw_do)/2;                  % center frequency
    bw=2*(bw_up-fc);
% 
%     figure(1);
%     plot(fx,pw,'k'); hold on;
%     plot(fx,fpw,'b');
%     plot([bw_up bw_up],[0 1],'g');
%     plot([bw_do bw_do],[0 1],'g');
%     plot([fc fc],[0 1],'k--');
%     title('Original signal');

    data=zeros(size(s.data));

    % band pass filtering
    disp('Band Pass filtering');
    filter_bounds=[max([0 fc-2*bw]) min([Fs/2 fc+4*bw])];
    filter_bounds=[filter_bounds(1) filter_bounds(1)+0.5e6 filter_bounds(2)-0.5e6 filter_bounds(2)];
    data = band_pass(s.data,Fs,filter_bounds);
% 
%     [fx pw] = power_spectrum(data,Fs);
%     figure(1);
%     plot(fx,pw,'r--'); 

    % demodulation
    dt=1/Fs;
    mod_sig=exp(-j*2*pi*fc*s.time)*ones(1,size(data,2)); % demodulation sognal
    wb = waitbar(0, 'Demodulating...');
    for f=1:size(data,4)
        for n=1:size(data,3)
            waitbar(n / size(data,3), wb);
            data(:,:,n,f)=data(:,:,n,f).*mod_sig;
        end
    end
    close(wb);
% 
%     [fx pw] = power_spectrum(data,Fs);
%     figure(2);
%     plot(fx,pw,'b-'); hold on;

    % low pass filtering
    disp('Base Band filtering');
    filter_bounds=[fc-0.5e6 fc];
    data=low_pass(data,Fs,filter_bounds);

%     [fx pw] = power_spectrum(data,Fs);
%     figure(2);
%     plot(fx,pw,'r--'); 

    % resampling
    temp=[];
    fs=round(4*bw/1e6)*1e6; % new sampling frequency
    dt=1/fs;
    t=(s.time(1):(1/fs):s.time(end));
    wb = waitbar(0, 'Resampling...');
    for f=1:size(data,4)
        for ntx=1:size(data,3)
            for nrx=1:size(data,2)
                temp(:,nrx,ntx,f)=interp1(s.time,data(:,nrx,ntx,f),t,'linear',0);
            end
            waitbar(ntx / size(data,3), wb);
        end
    end
    close(wb);
       

%     [fx pw] = power_spectrum(temp,fs);
%     figure(3);
%     plot(fx,pw,'k-'); 

    % write in the output
    out.modulation_frequency=fc;
    out.format=E.signal_format.IQ;
    out.time=t.';
    out.data=temp;
end








