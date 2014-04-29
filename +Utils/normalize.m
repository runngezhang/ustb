function x = normalize(y,method,dim)

if nargin == 1
    method = 'max';
end

if nargin == 2
    dim = 1;
end

D =  ones(1, max(ndims(y),dim));        
D(dim) = size(y,dim);
        
switch method
    case 'max'
        nf = repmat(max(y,[],dim),D);        
    case 'sum'                
        nf = repmat(sum(y,dim),D);                        
    otherwise
        nf = repmat(max(y,[],dim),D);
end

x = y./nf;