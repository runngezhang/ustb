function varargout=find_indx(x,lim)

ix = zeros(size(lim));
for k=1:length(lim)
    [dummy ix(k)] = min(abs(x-lim(k)));
end

if nargout >= 1 
    varargout{1} = ix;
end

if nargout == 2     
    varargout{2} = ix(1):ix(end);        
end



    