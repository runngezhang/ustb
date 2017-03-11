function out_dataset=complex_demodulation(h)
%% Hilbert tranform demodulation
% 2015-03-26 Øyvind K.-V. Standal
% 2015-06-04 Alfonso Rodriguez-Molares

% fill the gaps
[h out_dataset]=h.fill_gaps();

doresample = false; % need to interpolate (not just downsample by integer rate)
if isempty(h.downsample_frequency),
    Ndown = 1; % no downsampling
else
    Ndown = round(h.sampling_frequency/h.downsample_frequency);
    %if abs(Ndown - round(Ndown)) < 0.01,
    %    Ndown = round(Ndown);
    %else
    %    doresample = true;
    %    warning('Downsampling by a non-integer rate is highly inefficient.');
    %end
end

siz = size(h.raw_data.data);
N = prod(siz(3:end)); % collapse dimensions beyond 3
ind_new = 1:Ndown:siz(1);
L = length(ind_new); % length of downsampled array

% preallocate final array
iq = zeros([L,siz(2:end)], 'like', 1i);

% demodulation mixing vector
mix = exp(-2i*pi*h.modulation_frequency/h.sampling_frequency*(0:siz(1)-1)');

% otherwise not really need for low-pass filtering
dofilter = (Ndown > 1);
if dofilter, % make low-pass filter
    % calculate window parameters
    [M, Wn, beta, typ] = kaiserord([0.9,1]*h.modulation_frequency, [1,0], [1e-2,1e-2], h.sampling_frequency);
    b = fir1(M, Wn, typ, kaiser(M+1,beta), 'noscale'); % filter design
    assert(size(h.raw_data.data,1)>3*(length(b)-1),'Too many low-pass filter coefficients! Try relaxing the filter design, or increasing the samples in the dataset.');
end

% preallocate temp array for Hilbert
% H = zeros(siz(1:2), 'like', 1i);

for nn = 1:N, % loop over higher dimensions
    H = bsxfun(@times, hilbert(h.raw_data.data(:,:,nn)), mix);
    if dofilter,
        H(:,:) = filtfilt(b, 1, H); % lowpass filter
    end
    if doresample,
        iq(:,:,nn) = interp1(H, ind_new, 'linear', 0);
    elseif Ndown >= 2,
        iq(:,:,nn) = H(ind_new,:); % downsample
    else
        iq(:,:,nn) = H;
    end
end

% downsampling
dt=Ndown/h.sampling_frequency;
t=(h.raw_data.time(1):dt:h.raw_data.time(end));

% copy results
out_dataset.modulation_frequency=h.modulation_frequency;
out_dataset.time=t.';
out_dataset.raw_data.data=iq;
end