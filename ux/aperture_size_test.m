max_depth = 60; %mm
max_aperture_size = 46; %elements

aperture_curve = [15 53 100];
mid_point = 50;

focus_depth = 30; %mm

% compute HALF aperture size in number of elements
aperture_size_el = compute_aperture_size(aperture_curve,mid_point,max_depth,max_aperture_size,focus_depth);

% compute the complete aperture size. We subtract 1 instead of adding one
% since that is what ultrasonix do.
aperture_size_el = 2*aperture_size_el - 1;
