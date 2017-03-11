function out_dataset=original(h)
%% ORIGINAL: band-pass + demodulation + low pass filter
% 2015-05-15 Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>

% fill the gaps
[h out_dataset]=h.fill_gaps();

% power spectrum
[fx pw] = tools.power_spectrum(h.raw_data.data,h.sampling_frequency);
assert(sum(pw)>0,'Dataset is zero');

% computing central frequency and bandwidth
disp('Estimating power spectrum');
fpw=filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)];
[dc ic]=max(fpw.*(fx>0).'); fc=fx(ic);
bw_up=min(fx((fx>fc)&(fpw<dc/2).')); % -6dB upper limit
bw_do=max(fx((fx<fc)&(fpw<dc/2).')); % -6dB down limit
fc=(bw_up+bw_do)/2;                  % center frequency
bw=2*(bw_up-fc);

if(h.plot_on)
    assert(sum(pw)>0,'Dataset is zero');
    figure();
    subplot(1,2,1);
    plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
    plot([h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
    plot(-[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
    axis([-2*bw_up*1e-6 2*bw_up*1e-6 0 1]);
    xlabel('f [MHz]');
    ylabel('Relative amplitude');
    title('Before demodulation');
end

% band pass filtering
disp('Band Pass filtering');
if isempty(h.bandpass_frequency_vector)
    transition=bw/10;
    low_freq=max([0 fc-2*bw]);
    high_freq=min([h.sampling_frequency/2*0.99 fc+2*bw]);
    h.bandpass_frequency_vector=[low_freq low_freq+transition high_freq-transition high_freq];
end
[data fre_res w]= tools.band_pass(h.raw_data.data,h.sampling_frequency,h.bandpass_frequency_vector);

if(h.plot_on)
    [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
    plot(h.sampling_frequency*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
    plot(fx,pw,'r--');
    plot([h.modulation_frequency h.modulation_frequency],[0 1],'m:');
end

% demodulation
mod_sig=exp(-j*2*pi*h.modulation_frequency*h.raw_data.time)*ones(1,size(data,2)); % demodulation sognal
wb = waitbar(0, 'Demodulating...');
for f=1:size(data,4)
    for n=1:size(data,3)
        waitbar(n / size(data,3), wb);
        data(:,:,n,f)=data(:,:,n,f).*mod_sig;
    end
end
close(wb);

if(h.plot_on)
    [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
    subplot(1,2,2);
    plot(fx*1e-6,pw,'k-'); hold on; grid on; axis manual;
    axis([-2*bw_up*1e-6 2*bw_up*1e-6 0 1]);
    xlabel('f [MHz]');
    ylabel('Relative amplitude');
    title('After demodulation');
end

% low pass filtering
disp('Base Band filtering');
baseband_frequency_vector=[0.9*h.modulation_frequency h.modulation_frequency];
[data fre_res w]=tools.low_pass(data,h.sampling_frequency,baseband_frequency_vector);

if(h.plot_on)
    [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
    subplot(1,2,2);
    plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
    plot(h.sampling_frequency*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
    plot(-h.sampling_frequency*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
end

% resampling
dt=1/h.downsample_frequency;
t=(h.raw_data.time(1):dt:h.raw_data.time(end));
data=reshape(data,size(h.raw_data.data)); % we get back singleton dimmensions
downsample_data=zeros(length(t),size(data,2),size(data,3),size(data,4));
wb = waitbar(0, 'Resampling...');
for f=1:size(data,4)
    for ntx=1:size(data,3)
        for nrx=1:size(data,2)
            downsample_data(:,nrx,ntx,f)=interp1(h.raw_data.time,data(:,nrx,ntx,f),t,'linear',0);
        end
        waitbar(ntx / size(data,3), wb);
    end
end
close(wb);

if(h.plot_on)
    [fx pw] = tools.power_spectrum(downsample_data,h.downsample_frequency);
    plot(fx*1e-6,pw,'m:');
end

% write in the output
out_dataset.modulation_frequency=h.modulation_frequency;
out_dataset.initial_time=t(1);
out_dataset.sampling_frequency=1/(t(2)-t(1));
out_dataset.data=downsample_data;
end

