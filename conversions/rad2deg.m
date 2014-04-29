function y=rad2deg(x)

y = 360/2/pi;
if nargin == 1
    y = y*x;
end