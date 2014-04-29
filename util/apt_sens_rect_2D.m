function [He,PHI,PSI] = apt_sens_rect_2D(a,b,lambda,phi,psi)
% He = apt_sens_rect_2D(a,b,lambda,phi,psi)
%
% Computes the directivity sensitivity function for a rectangluar aperture
% according to equation 5.190 in Bjørns book.
%
% Input:
% a         - Azimuth aperture half size
% b         - Azimuth aperture half size
% lambda    - Wavelength
% phi       - Directions azimuth in radians
% psi       - Directions elevation in radians
%
% Output: 
% He        - Sensitivity function
% PHI       - 2D axis for plotting
% PSI       - 2D axis for plotting

[PHI PSI] = meshgrid(phi,psi);

He = sinc((2*a/lambda)*sin(PHI)).*sinc((2*b/lambda)*sin(PSI)).*cos(PHI).*cos(PSI);