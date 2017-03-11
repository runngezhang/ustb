function out_dataset=fast(h)
%% FAST: Demodulation + low pass filter
% 2015-09-14 Alfonso Rodriguez-Molares

% fill the gaps
[h out_dataset]=h.fill_gaps();

siz = size(h.raw_data.data);
N = prod(siz(3:end));        % collapse dimensions beyond 3

% copy data
data = h.raw_data.data(:,:,1:N);

% show spectrum
if(h.plot_on)
    [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
    assert(sum(pw)>0,'Dataset is zero');
    figure();
    subplot(1,2,1);
    plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
    plot([h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
    plot(-[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
    axis([-4*h.modulation_frequency*1e-6 4*h.modulation_frequency*1e-6 0 1]);
    xlabel('f [MHz]');
    ylabel('Relative amplitude');
    title('Before demodulation');
end

% demodulation
mod_sig=repmat(exp(-j*2*pi*h.modulation_frequency*h.raw_data.time),1,size(data,2),size(data,3));
data=data.*mod_sig;

if(h.plot_on)
    [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
    subplot(1,2,2);
    plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
    plot([0 0]*1e-6,[0 1],'r--');
    plot(-2*[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
    axis([-4*h.modulation_frequency*1e-6 4*h.modulation_frequency*1e-6 0 1]);
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
Ndown=round(h.sampling_frequency/h.downsample_frequency);
ind_new = 1:Ndown:siz(1);       % decimating vector
L = length(ind_new);            % length decimating  vector

% copy results
out_dataset.initial_time=h.raw_data.time(ind_new(1));
out_dataset.sampling_frequency=1/(h.raw_data.time(ind_new(2))-h.raw_data.time(ind_new(1)));
out_dataset.data=reshape(data(ind_new,:,:),[L siz(2:end)]);   % resampled data
out_dataset.modulation_frequency=h.modulation_frequency;
end
