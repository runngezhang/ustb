function [H r] = spat_freq_resp_annular(F,a,f,w,c0,rho,retarded_time)
% H = spat_freq_resp_annular(F,a,f,w,c0,rho)
%
% Computes the spatial frequency response for an annular array according to
% equation 5.254 in Volume I. The computed spatial frequency response has
% dimensions size(H) = [length(w) length(f)].
%
% Input:
% F     - Focus depth
% a     - Aperture radius
% f     - Frequency axis
% w     - spatial distance from axis in the focal plane
% c0    - Speed of sound
% rho   - Mass density
%
% Output:
% H     - Spatial frequency response
%
% Author: Thor Andreas Tangen, ITK, ISB.
% Last change: 04 / 07 / 2010.

Nfft = length(f);
Nw = length(w);
k = 2*pi*f/c0;

k = k(:,ones(1,Nw));
w = w(:,ones(1,Nfft))';    

r = sqrt(F*F + w(1,:)'.^2);    
r = r(:,ones(1,Nfft))';    

if retarded_time
    rt = F;
else
    rt = 0;
end

A = a*a*exp(-1i*k.*(r-rt))./(2*r);
    
D = (a*k.*w/F);
J1 = besselj(1,D);

H = rho*2*A.*(J1./D);
H(k == 0) = 0;

ss = find(w(1,:) == 0,1);

if ~isempty(ss);
    ix0 = find_indx(ss,1);
    H(:,ss) = rho*a*a*exp(-1i*k(:,1)*(F-rt))./(2*F);
end