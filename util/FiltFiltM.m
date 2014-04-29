function y = filtfilt(b,a,x)
%FILTFILT Zero-phase forward and reverse digital filtering.
%   Y = FILTFILT(B, A, X) filters the data in vector X with the filter described
%   by vectors A and B to create the filtered data Y.  The filter is described 
%   by the difference equation:
%
%     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%
%   After filtering in the forward direction, the filtered sequence is then 
%   reversed and run back through the filter; Y is the time reverse of the 
%   output of the second filtering operation.  The result has precisely zero 
%   phase distortion and magnitude modified by the square of the filter's 
%   magnitude response.  Care is taken to minimize startup and ending 
%   transients by matching initial conditions.
%
%   The length of the input x must be more than three times
%   the filter order, defined as max(length(b)-1,length(a)-1).
%
%   Note that FILTFILT should not be used with differentiator and Hilbert FIR
%   filters, since the operation of these filters depends heavily on their
%   phase response.
%
%   See also FILTER.

%   References: 
%     [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed., McGraw-Hill, 2001
%     [2] Fredrik Gustafsson, Determining the initial states in forward-backward 
%         filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
%         Volume 44, Issue 4

%   Author(s): L. Shure, 5-17-88
%   revised by T. Krauss, 1-21-94
%   Initial Conditions: Fredrik Gustafsson
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.4 $  $Date: 2009/08/11 15:47:37 $

    error(nargchk(3,3,nargin,'struct'))

% Jan Simon, 06-Aug-2009: Vectorized call of FILTER.
% The original source is property of The MathWorks!
% Patched: 23-Sep-2010 14:03:54

    if (isempty(b) || isempty(a) || isempty(x))
        y = [];
        return
    end

    m = size(x, 1);
    
    if m==1
        x = x(:);   % convert row to column
    end
    ndimsx = ndims(x);
    if ndimsx > 2
        saveSize = size(x);
        x = reshape(x,m,[]);
    end;
    len = size(x,1);   % length of input
    b = b(:).';
    a = a(:).';
    nb = length(b);
    na = length(a);
    nfilt = max(nb,na);

    nfact = 3*(nfilt-1);  % length of edge transients

    if (len<=nfact),    % input data too short!
        error(generatemsgid('InvalidDimensions'),'Data must have length more than 3 times filter order.');
    end

% set up filter's initial conditions to remove dc offset problems at the 
% beginning and end of the sequence
    if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
    if na < nfilt, a(nfilt)=0; end
% use sparse matrix to solve system of linear equations for initial conditions
% zi are the steady-state states of the filter b(z)/a(z) in the state-space 
% implementation of the 'filter' command.
    rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
    cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
    data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
    sp = sparse(rows,cols,data);
    zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
% non-sparse:
% zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
%      ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

% Extrapolate beginning and end of data sequence using a "reflection
% method".  Slopes of original and extrapolated sequences match at
% the end points.
% This reduces end effects.
lenOnes    = zeros(1, nfact);
lenOnes(:) = len;
y = [ 2 * x(ones(1, nfact), :) - x((nfact + 1):-1:2, :); ...
      x; ...
      2 * x(lenOnes,        :) - x((len - 1):-1:(len - nfact), :)];

% Filter, reverse data, filter again, and reverse data again:
y = filter(b, a, y, zi * y(1, :));
y = y(size(y, 1):-1:1, :);
y = filter(b, a, y, zi * y(1, :));

% Remove extrapolated pieces of y:
y = y((len + nfact):-1:(nfact + 1), :);
    if ndimsx>2
        y = reshape(y,saveSize);
    end;
    if m == 1
        y = y.';   % convert back to row if necessary
    end
   
