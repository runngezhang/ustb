function rf_f = rf_depth_filtering(rf,Ts,Ntaps,fco1,ix1,fco2,ix2,fco3,overlap)
% rf_f = rf_depth_filtering(rf,h1,ix1,h2,ix2,h3,overlap)
%
% Function for performing depth dependent RF filtering. 3 FIr filters are
% applied for 3 different depth regions
%
% Input:
% rf        - RF signal to be filtered
% Ts        - Sampling interval
%           - Number of filter coeffs 
% h1        - First RF filter cutoffs (Hz)
% ix1       - Depth limit for first RF filter (samples)
% h2        - Second RF filter cutoffs (Hz)
% ix2       - Depth limit for second RF filter (samples)
% h3        - Third RF filter cutoffs (Hz)
% ocverlap  - Overlap between depth regions (samples)
%
% Output:
% rf_f      - Filtered RF signal
% 
% Author : Thor Andreas Tangen
% Last edited : 2010-09-22 TAT

h1 = fir1(Ntaps,2*fco1*Ts,hamming(Ntaps+1));
h2 = fir1(Ntaps,2*fco2*Ts,hamming(Ntaps+1));
h3 = fir1(Ntaps,2*fco3*Ts,hamming(Ntaps+1));

dims = size(rf);

% Make 2-dimensional, easier to perform filtering
rf = rf(:,:);
rf_f = zeros(size(rf));

% Make overlap an even number and divide by 2
overlap = (overlap + mod(overlap,2))/2;
ix3 = size(rf,1);

weights = zeros(ix3,1);

% FIR 1
indxs = 1:(ix1+overlap);
rf_f(indxs,:) = filter(h1,1,rf(indxs,:),[],1);
weights(indxs) = weights(indxs) + 1;

% FIR 2
% Start length(h2) - 1 samples before in order for stabilizing FIR filter
indxs = (ix1-overlap-length(h2)+1):(ix2+overlap);
tmp = filter(h2,1,rf(indxs,:),[],1);
% Remove the length(h2) - 1 first samples
indxs = (ix1-overlap):(ix2+overlap);
rf_f(indxs,:) = rf_f(indxs,:) + tmp((length(h2)):end,:);
weights(indxs) = weights(indxs) + 1;

% FIR 3
% Start length(h3) - 1 samples before in order for stabilizing FIR filter
indxs = (ix2-overlap-length(h3)+1):ix3;
tmp = filter(h3,1,rf(indxs,:),[],1);
% Remove the length(h3) - 1 first samples
indxs = (ix2-overlap):ix3;
rf_f(indxs,:) = rf_f(indxs,:) + tmp((length(h3)):end,:);
weights(indxs) = weights(indxs) + 1;

weights(weights == 0) = 1;
rf_f = rf_f./weights(:,ones(1,size(rf_f,2)));

rf_f = reshape(rf_f,dims);
