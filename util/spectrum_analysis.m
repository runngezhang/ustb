function [ix fc B_dB rB_dB flim]=spectrum_analysis(G,f,lim,unit)
%SPECTRUM ANALYSIS
%[dB_ix fc B_6dB rB_6dB]=spectrum_analysis(G,f,lim,unit)
%Finds the indx of the given dB limit. If no dB limit is given 6dB is used.
%Input:
%G -  Energy spectrum, G = |fft(g)|
%f - frequency axis for the given spectrum
%lim - Optional, the limit to be used, specified in the units, default 6dB
%unit - Optional, unit of the limits, 'dB' ot 'Hz'
%
%Output:
% ix - index for the given limits
% fc - center frequency, computed using the dB interval specified, fc =
% \int_{f_{dB}}df fG(f)
% B_dB - the bandwidth in Hz, with the given db limits, only computed when
% unit is dB
% rB_dB - relative bandwidth, only computed when
% unit is dB
% flim - gives at what frequency the limits are. Interpolated, not the same as f(ix).

if nargin == 2
    lim = -6;
    unit = 'dB';   
end

if nargin == 3
    unit = 'dB';
end

if ~isreal(G)
    G = abs(G);
end

if strcmp(unit,'dB')
    [dB_ix flim] = find_dB_limits(G,lim);

    df = f(2) - f(1);

    P = sum(G(dB_ix(1):dB_ix(2))*df);
    
    if size(G(dB_ix(1):dB_ix(2))) ~= size(f(dB_ix(1):dB_ix(2)))
        f = f';
    end
    
    fc = (1/P)*sum(G(dB_ix(1):dB_ix(2)).*f(dB_ix(1):dB_ix(2))*df);
    %B_dB = f(dB_ix(2))-f(dB_ix(1));
    flim = interp1(1:length(f),f,flim,'linear');
    B_dB = flim(2) - flim(1);

    rB_dB = B_dB/fc;

    ix = dB_ix;
else
    df = f(2) - f(1);
    f_ix = find_indx(f,lim);
    P = sum(G(f_ix(1):f_ix(2))*df);
    
    if size(G(f_ix(1):f_ix(2))) ~= size(f(f_ix(1):f_ix(2)))
        f = f';
    end
    
    fc = (1/P)*sum(G(f_ix(1):f_ix(2)).*f(f_ix(1):f_ix(2))*df);
    ix =f_ix;
end

