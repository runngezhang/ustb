function [ti dt dz] = compute_focusing_delays_phased(r,theta,ac,c0)
% COMPUTE_FOCUSING_DELAYS_PHASED
% Function for computing focusing delays for phased arrays.
% [ti dt dz] = compute_focusing_delays_phased(r,theta,ac,c0)
% Input:
% r - focusing points in the direction of theta
% theta - steering angle
% ac - distance from the center element to the center of each of the other
% elements
% c0 - speed of sound

theta = 2*pi*theta/360;
xa = ac;

zF = r*cos(theta);
xF = r*sin(theta);

xF = xF(ones(1,length(xa)),:);
zF = zF(ones(1,length(xa)),:);
xa = xa(ones(1,length(r)),:)';

ri = sqrt((xF-xa).^2 + zF.^2);
ti = (r(ones(1,length(ac)),:) + ri)/c0;

if nargin > 1    
    dt = ti - 2*r(ones(1,length(ac)),:)/c0;
end

if nargout > 2
    dz = c0*dt/2;
end

