function [ti dt dz] = compute_focusing_delays_linear(z,ac,c0)

d1 = max(size(z));
d2 = max(size(ac));
z = z(ones(1,d2),:);
ac = ac(ones(1,d1),:).';

Rz = sqrt(ac.^2 + z.^2);
ti = (z + Rz)/c0;

if nargout > 1
    dt = (z - Rz)/c0;% ti - 2*z/c0;
end

if nargout > 2
    dz = c0*dt/2;
end