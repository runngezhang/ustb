function [dB_ix dB_x] = find_db_limits(G,lim,x,tight)
% function [dB_ix val] = find_db_limits(G,lim)
% function [dB_ix val] = find_db_limits(G,lim,x)
% function [dB_ix val] = find_db_limits(G,lim,x,tight)
% function [dB_ix val] = find_db_limits(G,lim,[],tight)
%
% Finds the db limits of the given data, G, relative to the maximum 
% value of G.
%
% The default behaviour is this:
% Starts searching from the first sample of G and finds the first data
% point above the specified db limit. From this sample it starts searching
% for the last data point which is above the db limit. 
%
% If the <i>tight</i> parameter is set to true, the behaviour is the
% following:
% Starts searching from the first sample of G and finds the first data
% point above the specified db limit. From this sample it starts searching
% for the <i>first</i> data point which is below the db limit. 
% 
% <i>dB_ix</i> is the indices found, while <i>dB_x</i> is the linearly
% interpolated points of the db limits.
%
% Input:
% G         - Data series. If imaginary, abs(G) is used.
% lim       - db limit, should be negative
% x         - (Optional) Sample points for the data series [default: 1:length(G)]
% tight     - (Optional) Specify the search behaviour, see above. [default: 0]
%
% Output:
% dB_ix     - Indices of the limit points
% dB_x      - Linearly interpolated points of the dB limits accoriding to
% the sample vector x
%
% Example:
% Nfft = 1024;
% X = abs(fft(x,Nfft));
% f = fs*(0:(Nfft-1))/Nfft;
% lim = -6;
% 
% [db_ix db_lim_f] = find_db_limits(G(1:Nfft/2),lim,f(1:Nfft/2));
%
% ylims = [-40 0];
% figure;
% plot(f,20*log10(G),...
%   db_lim_f([1 1]),ylims,'k--',...
%   db_lim_f([2 2]),ylims,'k--');
% xlim(f([1 Nfft/2])),ylim(ylims)
% grid on
%
% Author: Thor Andreas Tangen, InPhase Solutions AS
% Copyright InPhase Solutions AS
% Created 20.10.2010
% 

    G = shiftdim(G);
    
    if nargin < 4
        tight = 0;
    end
    
    if nargin < 3 || isempty(x)
        x = 1:size(G,1);
    end
    
    if lim > 0
        lim = -lim;
    end

    if ~isreal(G)
        G = abs(G);
    end
    
    [R C] = size(G);
    dB_ix = zeros(C,2);
    dB_x = [0 0];
    for c=1:C
        G_dB = 20*log10(G(:,c)./max(G(:,c)));

        [~, max_ix] = max(G_dB);

        ix = find((G_dB(1:max_ix)-lim) >= 0,1,'first'); 
        if isempty(ix)
            ix = 1;
        end
        dB_ix(c,1) = ix;
        if tight
            ix = find((G_dB((max_ix+1):end)-lim) <= 0,1,'first');
        else
            ix = find((G_dB((max_ix+1):end)-lim) >= 0,1,'last');
        end
        if isempty(ix)
            ix = R-max_ix+1;
        end

        dB_ix(c,2) = min(ix + max_ix+1,length(x));        
        
        if nargout == 2

            if nargin < 3
                x = (0:(length(G)-1))';
            end
            if size(x,2) == 1
                x = x';
            end
            try               
                if dB_ix(c,1) == 1
                    dB_x(c,1) = 1;
                else
                    a = G_dB(dB_ix(c,1))-G_dB(dB_ix(c,1)-1);
                    b = G_dB(dB_ix(c,1)-1);
                    xx = (lim-b)./a;
                    dB_x(c,1) = x(dB_ix(c,1)-1) + (x(dB_ix(c,1))-x(dB_ix(c,1)-1)).*xx';
                end
                
                a = G_dB(dB_ix(c,2))-G_dB(dB_ix(c,2)-1);
                b = G_dB(dB_ix(c,2)-1);
                if a == 0
                    dB_x(c,2) = x(dB_ix(c,2));
                else
                    xx = (lim-b)./a;
                    dB_x(c,2) = x(dB_ix(c,2)-1) + (x(dB_ix(c,2))-x(dB_ix(c,2)-1)).*xx';
                end
            catch e
                rethrow(e);
            end
        end
    end
end