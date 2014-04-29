function [xi,yi, d2x] = interp_max(y,xs,method)
% [xi,yi] = interp_max(y,xs,method)
%
% Function for for finding the subsample maximum of a dataseries by
% interpolation. The function uses parabolic interpolation around the max
% of the dataseries to fin the subsample maximum value and origin of the
% maximum value.
%
% Input:
% x         - Data series
% xs        - Sampling interval. Optional, default xs = 1
% method    - Method, 'parabolic' or 'cubic'. Optional, 'parabolic' default
%
% Output:
% xi        - Origin of max value
% yi        - Max value
if nargin < 2 || isempty(xs)
    xs = 1;
end    

if nargin < 3
    method = 'parabolic';
end

[dummy ix] = max(y);

switch method
    case 'parabolic'
        if ix == 1 || ix == length(y)
            xi = xs*ix;
            return
        end
        
        xi = xs*((y(ix-1) - y(ix+1))/(2*(y(ix-1)-2*y(ix)+y(ix+1))) + ix - 1);   
                
        if xi > xs*(ix-2) && xi < xs*ix
            six = ix-1;
        else
            six = ix;
        end
        
        c = y(six);
        b = (y(six+1) - y(six))/xs;
        a = ((y(six+2) - y(six+1))/xs - (y(six+1) - y(six))/xs)/(2*xs);        
        yi = c + b*(xi - xs*six) + a*(xi - xs*six)*(xi - xs*(six +1));        
        d2x = 2*a;        
    case 'cubic'
        ixs = (-1:2)';
        if size(y,1) == 1
            y = y';
        end
        pp = polyfit(ixs,y(ix + ixs),3);
        
        a = 3*pp(1);
        b = 2*pp(2);
        c = pp(3);
        xi = (-b - sqrt(b*b - 4*a*c))/2/a;
        
        yi = pp(1)*xi*xi*xi + pp(2)*xi*xi + pp(3)*xi + pp(4);
        xi = xs*(ix + xi - 1);                
        d2x = 6*p(1)*xi + 2*p(2);
end
 
