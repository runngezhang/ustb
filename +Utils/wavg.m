function ex = wavg(x,w,dx,dim)
% EX = WAVG(X,W)
% Finds the weighted average of X. W is the weigths. If X is a matrix the
% weighted average is done columnwise. W can then be either the same size
% as X or length(W) = size(X,1).

if nargin < 4
    dim = 1;
end

if nargin < 3
    dx = 1;
end


%w = w./sum(w*dx);
ex = sum(x.*w*dx,dim)./sum(w*dx,dim);

