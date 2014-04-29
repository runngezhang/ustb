 function [dB_ix val] = find_dB_limits(G,lim,x)

    if lim > 0
        lim = -lim;
    end

    if ~isreal(G)
        G = abs(G);
    end

    G = shiftdim(G);

    [R C] = size(G);
    dB_ix = zeros(C,2);
    val = [0 0];
    for c=1:C
        G_dB = 20*log10(G(:,c)./max(G(:,c)));

        [dummy max_ix] = max(G_dB);

        ix = max(find((G_dB(1:max_ix)-lim) > 0,1,'first'),2); 
        if isempty(ix)
            ix = 1;
        end
        dB_ix(c,1) = ix;

        ix = find((G_dB(max_ix:end)-lim) > 0,1,'last');
        if isempty(ix)
            ix = R-max_ix+1;
        end

        dB_ix(c,2) = min(ix + max_ix,length(x));        
        if nargout == 2

            if nargin < 3
                x = (0:(length(G)-1))';
            end
            if size(x,2) == 1
                x = x';
            end
            try
                %dx = x(2) - x(1);
                a = G_dB(dB_ix(c,:))-G_dB(dB_ix(c,:)-1);
                b = G_dB(dB_ix(c,:)-1);
                xx = (lim-b)./a;
                val(c,:) = x(dB_ix(c,:)-1) + (x(dB_ix(c,:))-x(dB_ix(c,:)-1)).*xx';
            catch e
                keyboard
            end
        end
    end
end