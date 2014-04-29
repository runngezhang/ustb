function logP=imagelog2(P,gain,dyn,handles)

%IMAGELOG   Image a matrix of ultrasound power in log scale
%
%   function logP=imagelog(P,gain,dyn);
%
%   Inputs: P    - Power image
%           gain - gain [dB] (=image brightness)
%           dyn  - dynamic range [dB] (=image contrast)
%
%   Output: logP - Output image, log Power.
%
%   If no output image is requested, the image is displayed
%   in 64 graylevels in the current axes. If the input image is
%   complex valued, only the real part is displayed.
%
% NB! image values should be from 1 to length(colormap)
% NOT 0 to length(colormap)-1! See line 47

% if nargin<3, error('Too few input arguments');end;
% 
% if isstruct(handles),
%     maxIndex = handles.cmapmax-1;
% else
%     maxIndex = handles;
% end
% 
% if ~isreal(P),
%     disp('Warning: Imaginary part of image discarded');
% end;
% 
% if prod(size(gain))~=1|~isreal(gain)|isstr(gain),
%     error('Invalid gain. Must be a real number');
% end;
% 
% if prod(size(dyn))~=1|~isreal(dyn)|dyn<=0|isstr(dyn),
%     error('Invalid dynamic range. Must be real and positive');
% end;

minIndex=0;
if nargin > 3 && isnumeric(handles)
    maxIndex = handles;
else
    maxIndex=255;
end
SmallNonZeroNumber=1e-10;%10^(-maxIndex);
% logP=maxIndex*(gain+10*log10(max(SmallNonZeroNumber,real(P))))/dyn;
% logP=maxIndex*(gain+10*log10(1+P))/dyn;
logP=maxIndex*gain/dyn+maxIndex/dyn*10*log10(max(SmallNonZeroNumber,real(P)));
logP = round(min(maxIndex, max(minIndex,logP)));
if ~nargout,
    image(logP);
    colormap(gray(maxIndex));
    logP=[];
end;