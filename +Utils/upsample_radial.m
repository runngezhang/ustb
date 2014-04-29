function x_us = upsample_radial(x,N)
% x_us = upsample_radial(x)
%
% Function for doing upsampling by a factor of N in the radial direction.
% Should be used with IQ data in get best result. Uses cubic
% interpolation: 
% 
%       a0 = y3 - y2 - y0 + y1;
%       a1 = y0 - y1 - a0;
%       a2 = y2 - y0;
%       a3 = y1;
%       
%       ti = 1/N;
%       yi = a0*ti^3+a1*ti^2+a2**ti+a3;
%
% Input:
% x     - Input data
% N     - Order of upsampling
%
% Output:
% x_us  - Upsampled data

x_us = zeros(N*size(x,1),size(x,2),size(x,3));
x_us(1:N:end,:,:) = x(:,:,:);

x3 = x(4:end,:,:);
x2 = x(3:end-1,:,:);
x1 = x(2:end-2,:,:);
x0 = x(1:end-3,:,:);

a0 = x3-x2-x0+x1;
a1 = x0-x1-a0;
a2 = x2-x0;
a3 = x1;

t = 1/N;
for kk=1:(N-1)
    mu = kk*t;
    xi = ((a0*mu + a1)*mu + a2)*mu + a3;
    
    ix = (4+N+kk-3):N:size(x_us,2);
    x_us(ix(1:size(xi,1)),:,:) = xi;
end