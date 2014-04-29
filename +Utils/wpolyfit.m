function p = wpolyfit(x,y,n,w)
% WPOLYFIT
% Performs weighted polynomial fitting. Output can be used with other
% polynomial functions in matlab like polyval etc...
% TODO: SVD implementation of pseudo inverse, or other more numerical
% stable solution of the normal equations.
% Input
% x - time steps
% y - data
% n - polynomial order
% w - (optional) weights

if size(x,1) == 1
    x = x';
end

if size(x) ~= size(y)
    y = y';
end

if nargin == 3
    w = ones(size(x));
end

X = zeros(length(x),n+1);
for k=n:-1:0
    X(:,n-k+1) = x.^k;
end

W = diag(w);

p = inv(X'*W*X)*X'*W*y;
p = p';

