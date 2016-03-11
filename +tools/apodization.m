function apod = apodization(distance,aperture,window)
% function apodization = apodization(distance,aperture,window)
%
%   Assigns different apodization to a set of pixels and elements
%
%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2015/02/15 $

    switch(window)
        case 'none' 
            apod = ones(size(distance)); 
        case 'boxcar' 
            apod = double(distance<=aperture/2); 
        case 'hanning'
            apod = double(distance<=aperture/2).*(0.5 + 0.5*cos(2*pi*distance./aperture)); 
        case 'hamming'
            apod = double(distance<=aperture/2).*(0.53836 + 0.46164*cos(2*pi*distance./aperture)); 
        case 'tukey25'
            roll=0.25;
            apod =(distance<(aperture/2*(1-roll))) + (distance>(aperture/2*(1-roll))).*(distance<(aperture/2)).*0.5.*(1+cos(2*pi/roll*(distance./aperture-roll/2-1/2)));                               
        case 'tukey50'
            roll=0.5;
            apod=(distance<(aperture/2*(1-roll))) + (distance>(aperture/2*(1-roll))).*(distance<(aperture/2)).*0.5.*(1+cos(2*pi/roll*(distance./aperture-roll/2-1/2)));                               
        case 'tukey75'
            roll=0.75;
            apod=(distance<(aperture/2*(1-roll))) + (distance>(aperture/2*(1-roll))).*(distance<(aperture/2)).*0.5.*(1+cos(2*pi/roll*(distance./aperture-roll/2-1/2)));                               
        otherwise
            error('Unknown window type. Known types are: boxcar, hamming, hanning, tukey25, tukey50, tukey75.');
        end
end