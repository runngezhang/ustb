function [h X] = imagesc_logdyn(varargin)


if nargin == 1,% Only data given
   min_dyn = inf;
   dx = 1:size(varargin{1},2);
   dy = 1:size(varargin{1},1);
   X = varargin{1};
end
if nargin == 2,% Data and dynamic range given
   dx = 1:size(varargin{1},2);
   dy = 1:size(varargin{1},1);
   X = varargin{1};
   min_dyn = varargin{2};
end
if nargin == 3,%Dynamic range not given
   dx = varargin{1};
   dy = varargin{2};
   X = varargin{3};
   min_dyn = inf;
end
if nargin == 4,
   dx = varargin{1};
   dy = varargin{2};
   X = varargin{3};
   min_dyn = abs(varargin{4});
end

X = X/max(abs(X(:)));
X = 20*log10(max(abs(X),eps));
X(find(X<-min_dyn)) = -min_dyn;

h = imagesc(dx,dy,X);
colormap('Gray');
%axis image;
