function zc_ix = zerocrossings(x,type)
% ZEROCROSSINGS Finds the zero crossings of vector x
% For a sequence [6 5 3 -1 -2] the zero crossing index will be 3.5, while
% for a sequence [6 5 3 0 -1 -2], the zero crossing index will be 4. A
% sequence [6 5 3 0 2] will not produce a zero crossing at index 4. Platous
% like [6 5 3 0 0 0 -2] will give zerocrossings at
if nargin < 2
    type = 'all';
end

s = sign(x);
t = filter([1 1],1,s);
ix = 1;

zc_ix = -1*ones(size(x));
switch type
    case 'all'
        while ix < (length(x)-1)
            if t(ix)==0 && s(ix) ~= 0        
                zc_ix(ix) = ix-0.5;
                ix = ix + 1;
            elseif s(ix) == 0 && ix > 1
                ix1 = ix;
                at_end = false;
                while s(ix) == 0
                    ix = ix + 1;
                    if ix == length(x) &&  s(ix) == 0
                        at_end = true;
                        break
                    end
                end
                if ~at_end && (s(ix1-1)+s(ix)) == 0
                    zc_ix(ix) = (ix1 + ix-1)/2;
                end
            else
                ix = ix + 1;
            end 
        end
    case 'max'
        while ix < (length(x)-1)
            if t(ix)==0 && s(ix) ~= 0 && s(ix) < 0
                zc_ix(ix) = ix-0.5;
                ix = ix + 1;
            elseif s(ix) == 0 && ix > 1
                ix1 = ix;
                at_end = false;
                while s(ix) == 0
                    ix = ix + 1;
                    if ix == length(x) &&  s(ix) == 0
                        at_end = true;
                        break
                    end
                end
                if ~at_end && (s(ix1-1)+s(ix)) == 0
                    zc_ix(ix) = (ix1 + ix-1)/2;
                end
            else
                ix = ix + 1;
            end 
        end
end

zc_ix = zc_ix(zc_ix > 0);