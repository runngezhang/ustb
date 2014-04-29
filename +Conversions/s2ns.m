function x=s2ns(y)
    
x = 1e9;
if nargin == 1
    if ~iscell(y)
        x = x*y;
    else
        result = cell(size(y));
        for i = 1:numel(y)
            result{i} = y{i}*x;
        end;
        x = result;
    end;
end