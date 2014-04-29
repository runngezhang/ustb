function x_us = upsample(x,N,dim,method)
% x_us = upsample(x,N,dim)
%
% Function for doing upsampling by a factor of N in the dim dimension.
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
% N     - Order of upsampling. Default = 2.
% dim   - Dimension, should be 1 or 2. Default = 1.
% method    - 'cubic' or 'linear'
%
% Output:
% x_us  - Upsampled data
if nargin < 2
    N = 2;
end

if nargin < 3
    dim = 2;
end

if nargin < 4
    method = 'cubic';
end

dims = size(x);
if dim > length(dims)
    error('dim is greater then the number of dimensions in input data')
end

d = 1:length(dims);
d = [d(dim) d(d~=dim)];

x = permute(x,d);
x = x([ones(1,5),1:end,end*ones(1,5)],:);

dims_us = [N*dims(d(1)) dims(d(2:end))];
x_us = zeros(dims_us(1)+N*10,prod(dims_us(2:end)));

x_us(1:N:end,:) = x(:,:);

switch method
    case 'cubic'
        x3 = x(4:end,:);
        x2 = x(3:end-1,:);
        x1 = x(2:end-2,:);
        x0 = x(1:end-3,:);

        a0 = x3-x2-x0+x1;
        a1 = x0-x1-a0;
        a2 = x2-x0;
        a3 = x1;

        t = 1/N;
        for kk=1:(N-1)
            mu = kk*t;
            xi = ((a0*mu + a1)*mu + a2)*mu + a3;
            %xi = a0*mu*mu*mu + a1*mu*mu + a2*mu + a3;
            ix = (1+N+kk):N:size(x_us,1);
            x_us(ix(1:size(xi,1)),:) = xi;
        end
    case 'linear'        
        x0 = x(1:end-1,:);
        x1 = x(2:end,:);
        
        for kk=1:(N-1)            
            xi = x0 + (kk/N)*(x1 - x0);
            ix = (1+kk):N:(size(x_us,1));            
            x_us(ix(1:size(xi,1)),:) = xi;
        end    
    case 'linear_abs_angle'        
        x0 = abs(x(1:end-1,:));
        x1 = abs(x(2:end,:));
        x0_ang = angle(x(1:end-1,:));
        x1_ang = angle(x(2:end,:));
        
        for kk=1:(N-1)            
            xi = x0 + (kk/N)*(x1 - x0);
            xi_ang = x0_ang + (kk/N)*(x1_ang - x0_ang);
            xi = xi.*exp(1i*xi_ang);
            ix = (1+kk):N:(size(x_us,1));            
            x_us(ix(1:size(xi,1)),:) = xi;
        end    
end
x_us = x_us(N*5+1:(end-N*5),:);
x_us = reshape(x_us,dims_us);
x_us = ipermute(x_us,d);
