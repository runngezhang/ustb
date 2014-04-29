function [pos MaxT_HF MaxT_LF]= find_hf_position(p_hf,p_lf,lf_polarity,Ts,Tp_LF)
% [pos MaxT_HF MaxT_LF]= find_hf_position(hf_env,p_lf,lf_polarity,Ts,p_hf)
%
% Function for finding the position of the HF relative to the LF. Finds 
% the maximum of the hf. Then find the maximum of of lf_polarity*p_LF
% within +/-Tp_LF/2 of the max point of the HF. The phase relation is the 
% difference between the maximum points for the two pulses.
%
% Input:
% p_HF          - HF pulse or envelope
% p_lf          - The LF pulse
% lf_polarity   - The polarity og the LF pulse, 1 or -1
% Ts            - Sampling interval
% Tp_LF         - Period of LF pulse
%
% Output:
% pos       - The phase relation between the HF and LF in seconds
% MaxT_HF   - The time of maximum point of the HF
% MaxT_LF   - The time of maximum point of the LF

t = Ts*(0:(length(p_lf)-1));
MaxT_HF = interp_max(p_hf,Ts,'parabolic');
%MaxT_HF = sum(t'.*p_hf)./sum(p_hf);

[ix ixs] = find_indx(t,MaxT_HF + Tp_LF/2*[-1 1]); 

tmp = zeros(size(p_lf));
tmp(ixs) = p_lf(ixs);

MaxT_LF = interp_max(lf_polarity*tmp,Ts,'parabolic');
pos = -(MaxT_LF - MaxT_HF);
% figure(100);
% plot(t,p_hf,MaxT_HF*[1 1],[0 0.5]);
% pause;


