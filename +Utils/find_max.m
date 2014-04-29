function [max_ix max_x max_y d2x] = find_max(y,loc_max_lim)

dY = diff(y);
zc = floor(zerocrossings(dY,'max')) + 1;

max_ix = zeros(size(zc));
max_x = zeros(size(zc));
max_y = zeros(size(zc));
d2x = zeros(size(zc));

valid = true(size(zc));
for kk=1:length(zc)
    if zc(kk) == 1 || zc(kk) == length(y)
        max_ix(kk) = zc(kk);
    else
        if y(zc(kk)-1) > y(zc(kk))
            max_ix(kk) = zc(kk) - 1;
        elseif y(zc(kk)+1) > y(zc(kk))
            max_ix(kk) = zc(kk) + 1;
        else            
            max_ix(kk) = zc(kk);
        end
        
        ixs = max_ix(kk) + (-loc_max_lim:loc_max_lim);
        ixs = ixs(ixs >= 1 & ixs <= length(y));
        if any(y(ixs) > y(max_ix(kk)))
            valid(kk) = false;
        else           
%             try
                [max_x(kk) max_y(kk) d2x(kk)] = interp_max(y(ixs),1,'parabolic');            
                max_x(kk) = max_ix(kk) + (max_x(kk)-(length(ixs)+1)/2);            
%             catch e
%                 keyboard
%             end
        end
    end
end
max_ix = max_ix(valid);
max_x = max_x(valid);
max_y = max_y(valid);