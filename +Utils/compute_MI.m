function [MI_HF MI_LF MI_comb] = compute_MI(p_HF,p_LF,fc_HF,fc_LF)
% [MI_LF MI_HF MI_comb] = compute_MI(p_HF,p_LF,fc_HF,fc_LF)
%
% Computes MI for SURF pulse complexes. The combined MI is computed as
%
% MI_comb = (p_LF_min + p_LF_min)/sqrt(fc_HF)
% MI_comb = max(MI_comb , MI_LF , MI_HF);
%
% Input:
% p_HF  - HF pulse in MPa
% p_LF  - LF pulse in MPa
% fc_HF - Center frequency for the HF in Hz
% fc_LF - Center frequency for the LF in Hz
%
% Output:
% MI_HF - MI for HF pulse
% MI_LF - MI for LF pulse
% MI_comb - MI for combined pulse

p_LF_min = abs(min(p_LF));
p_HF_min = abs(min(p_HF));

MI_LF = p_LF_min/sqrt(fc_LF*1e-6);
MI_HF = p_HF_min/sqrt(fc_HF*1e-6);

MI_comb = (p_LF_min + p_HF_min)/sqrt(fc_HF*1e-6);

%MI_comb = max([MI_LF,MI_HF,MI_comb]);