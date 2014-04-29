function y=deg2rad(x)

y = 2*pi/360;
if nargin == 1
    y = y*x;
end