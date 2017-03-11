function out_dataset=cheap(h)
%% CHEAP 2-samples-per_wavelength demodulation
% authors: Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
% $2015/09/14$

% fill the gaps
[h out_dataset]=h.fill_gaps();

warning('The demodulation frequency in the "cheap" demodulation is h.sampling_frequency/4');
if ~isempty(h.modulation_frequency)
    assert(h.modulation_frequency==h.sampling_frequency/4,'The demodulation frequency in the "cheap" demodulation must be h.sampling_frequency/4');
end

siz = size(h.data);
N = prod(siz(3:end));        % collapse dimensions beyond 3

%% copy data
data = h.data(:,:,1:N);

%% show spectrum
if(h.plot_on)
    %fft_data = fftshift(fft(data(:,:), 1024));
    [fx pw] = tools.power_spectrum(h.data,h.sampling_frequency);
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

%% demodulation
L=siz(1)-mod(siz(1),4);
demod_data = zeros(size(data(1:4:L,:,:)));
demod_data(1:L/4,:,:) =     (data(1:4:L,:,:) + data(2:4:L,:,:)/j);
%demod_data(1:2:L/2,:,:) =     data(1:4:L,:,:) + exp(-j/5)*data(2:4:L,:,:)/j;
%demod_data(2:2:L/2,:,:) = - ( data(3:4:L,:,:) + exp(-j/5)*data(4:4:L,:,:)/j);

if(h.plot_on)
    [fx pw] = tools.power_spectrum(demod_data,h.sampling_frequency/2);
    subplot(1,2,2);
    plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
    plot([0 0]*1e-6,[0 1],'r--');
    plot(-2*[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
    axis([-4*h.modulation_frequency*1e-6 4*h.modulation_frequency*1e-6 0 1]);
    xlabel('f [MHz]');
    ylabel('Relative amplitude');
    title('After demodulation');
end

% copy dataset
out_dataset.time=h.time(1:4:L);                          % resampled time vector
out_dataset.data=reshape(demod_data,[L/4 siz(2:end)]);   % resampled data
out_dataset.modulation_frequency=h.sampling_frequency/4;
end
