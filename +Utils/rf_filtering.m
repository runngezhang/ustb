function [rf_f h] = rf_filtering(rf,Ntaps,fco,Ts,zero_phase)

if nargin < 5
    zero_phase = 1;
end
h    = fir1(Ntaps,2*fco*Ts,hamming(Ntaps+1));

if zero_phase
    dims = size(rf);
    rf_f = filtfilt(h,1,rf(:));
    rf_f = reshape(rf_f,dims);
else
    dims = size(rf);
    rf_f = filter(h,1,rf(:,:),[],1);
    %rf_f = FiltFiltM(h,1,rf(:));
    rf_f = reshape(rf_f,dims);
end