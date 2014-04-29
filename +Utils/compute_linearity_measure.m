function [tau_linearity,zi tau_diff] = compute_linearity_measure(z,tau,fc)
dz = 0.25;
zi = 0:dz:z(end);

tau_i = interp1(z,tau,zi,'spline');

tau_linearity = ones(floor(length(zi)/2),length(zi));
tau_diff = zeros(floor(length(zi)/2),length(zi));

wc = 2*pi*fc;
for kk=2:length(zi)        
    tau_th = interp1(zi,tau_i,zi(kk)/2,'linear');

    for ii=1:ceil(kk/2)
        tau1 = (tau_i(kk-ii) + tau_i(ii))/2;
        tau_linearity(ii,kk) = abs(sin(wc*(tau1 - tau_th)));
        tau_diff(ii,kk) = tau1 - tau_th;
    end
end