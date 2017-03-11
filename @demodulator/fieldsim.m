function out_dataset=fieldsim(h)
%% Fieldsim Hilbert tranform demodulation

% fill the gaps
[h out_dataset]=h.fill_gaps();

% demodulation
start_time = h.raw_data.time(1);
size_in    = size(h.raw_data.data);

% Find complex envelope (demodulate pre-env to zero)
data_out      = hilbert(h.raw_data.data(:,:)); % Find pre-envelope
demodVect     = exp(-1i*2*pi*h.modulation_frequency*(start_time + (0:size(data_out,1)-1)*1/h.sampling_frequency)).';
data_out      = reshape(bsxfun(@times, data_out, demodVect), size_in);

% downsampling
if isempty(h.downsample_frequency)
    warning('The downsampling frequency is not specified. Using 4*modulation_frequency');
    h.downsample_frequency=4*h.modulation_frequency;
end
t_new         = 0 : 1/h.downsample_frequency : h.time(end);
siz = size(data_out);
data_out      = interp1(h.time, data_out(:,:), t_new, 'linear',0);
data_out      = reshape(data_out, [size(data_out,1), siz(2:end)]);

% update dataset
out_dataset.time=t_new.';
out_dataset.data=data_out;
out_dataset.modulation_frequency=h.modulation_frequency;
end

