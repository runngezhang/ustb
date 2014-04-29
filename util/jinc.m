function y = jinc(x)
% y = jinc(x)
%
% Computes 2*J_1(pi*x)/pi*x, where J_1() is the bessel function of first
% kind. x == 0, produces 1, like sinc.
%
% Input:
% x     - Input
%
% Output:
% y     - Output

y = 2*besselj(1,pi*x)./(pi*x);
y(x == 0) = 1;

