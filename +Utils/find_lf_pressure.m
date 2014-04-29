function [exp_LF exp_LF_grad max_LF delta_exp_LF]= find_lf_pressure(env_hf,p_lf,Ts,method,Tp_HF)
% [exp_LF exp_LF_grad max_LF delta_exp_LF] = find_lf_pressure(env_hf,p_lf,Ts)
%
% Copmutes the epxerienced LF pressure, the max LF pressure, the
% experienced LF gradient and the largest LF pressure difference within the
% HF pulse. Uses weighted average over the hf envelope to estimate
% experienced HF pressure and gradient.
%
% Input:
% env_hf    - HF envelope
% p_lf      - LF pulse
% Ts        - Sampling interval
%
% Output:
% exp_LF        - Experienced LF pressure
% exp_LF_grad   - Experienced LF gradient
% max_LF        - Maximum LF pressure
% delta_exp_LF  - Maximum LF pressure difference within HF pulse

if nargin < 4
    method = 'max';
end
indx = find_db_limits(env_hf,-6);
indxs = indx(1):indx(2);

max_exp_lf = max(p_lf(indxs));
min_exp_lf = min(p_lf(indxs));
delta_exp_LF = max_exp_lf - min_exp_lf;

dp_lf = zeros(size(p_lf));
dp_lf(2:end-1) = (p_lf(3:end) - p_lf(1:end-2))/Ts/2;

switch method
    case 'max'
        [dummy, ix] = max(env_hf);
        exp_LF = p_lf(ix);%sum(env_hf.*p_lf)/sum(env_hf);
        exp_LF_grad = dp_lf(ix-1);
    case 'wavg'
        if nargin < 5
            ixs = 1:length(p_lf);
        else
            [dummy, ix] = max(env_hf);
            win_l = ceil(Tp_HF/Ts);
            ixs = ix + (-win_l:win_l);
        end
        
        exp_LF = sum(env_hf(ixs).*p_lf(ixs))/sum(env_hf(ixs));
        exp_LF_grad = sum(env_hf(ixs).*dp_lf(ixs))/sum(env_hf(ixs));
    case 'avg'
        if nargin < 5
            ixs = 1:length(p_lf);
        else
            [dummy, ix] = max(env_hf);
            win_l = ceil(Tp_HF/Ts);
            ixs = ix + (-win_l:win_l);
        end
        
        exp_LF = mean(p_lf(ixs));
        exp_LF_grad = mean(dp_lf(ixs));    
end


max_LF = max(p_lf);

% t = Ts*(0:(length(p_lf)-1));
% figure(10)
% plot(t,env_hf,t(indx([1 1])),0.3*[-1 1],t(indx([end end])),0.3*[-1 1],t,p_lf);