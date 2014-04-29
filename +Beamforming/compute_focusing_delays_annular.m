function [ti dt dz] = compute_focusing_delays_annular(z0,ac,c0,F)

if F == Inf
    e = zeros(size(ac));
else
    e = F - sqrt(F^2 - ac.^2);
end

z0 = z0(ones(1,length(e)),:);
e = e(ones(1,length(z0)),:)';
ac = ac(ones(1,length(z0)),:)';

Rz = sqrt((z0 - e).^2 + ac.^2);
ti = (z0 + Rz)/c0;

if nargout > 1
    dt = ti - 2*z0/c0;
end

if nargout > 2
    dz = c0*dt/2;
end





