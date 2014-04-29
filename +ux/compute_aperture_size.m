function aperture_size=compute_aperture_size(aperture_curve,mid_point,max_depth,max_aperture_size,focus_depth)
% COMPUTE_APERTURE_SIZE(aperture_curve,mid_point,max_depth,max_aperture_siz
% e,focus_depth)
%
% Parameters:
% aperture_curve        y (size axis) value of aperture curve. 3 point curve,
% given in percent
% mid_point             x (depth axis) value of mid point of curve in percent. 1st and 3rd point is 0 and
% 100 percent.
% max_depth             What depth 100 percent in x direction corresponds to
% max_aperture_size     What aperture size 100 percent in y direction
% corresponds to
% focus_depth           What depth to compute the aperture size for. Given
% in same dimension as max_depth

p1 = [0 curve_pos(aperture_curve(1),max_aperture_size)];
p2 = [curve_pos(mid_point,max_depth) curve_pos(aperture_curve(2),max_aperture_size)];
p3 = [max_depth curve_pos(aperture_curve(3),max_aperture_size)];

F = focus_depth;

slope = zeros(size(F));
y_intercept = zeros(size(F));

slope(p2(1) > F) = (p2(2) - p1(2))/(p2(1) - p1(1));
slope(p2(1) <= F) = (p3(2) - p2(2))/(p3(1) - p2(1));

y_intercept(p2(1) > F) = -p1(2);
y_intercept(p2(1) <= F) = slope(p2(1) <= F)*p2(1) - p2(2);

aperture_size = min(slope.*F - y_intercept,max_aperture_size);


function xi=curve_pos(x,x_max)

frac = min(x/100.0,1.0);	
xi = frac*x_max;