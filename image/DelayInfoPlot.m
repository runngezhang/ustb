function [ax, figh] = DelayInfoPlot(aDelayData, aDelayBounds, aMeanXBounds, figh)
    if ~exist('aDelayBounds','var')
        aDelayBounds = [];
    end;
    if ~exist('aMeanXBounds','var')
        aMeanXBounds = [];
    end;
    if ~exist('figh','var')
        figh = figure;
    else
        figure(figh);
    end;
            
    % show estimated delay image
    ax(1) = subplot(221);
    if ~isempty(aDelayBounds)
        imagesc(s2ns(aDelayData),s2ns(aDelayBounds));
    else
        imagesc(s2ns(aDelayData));
    end;        
    title('delay image');
    xlabel('x');
    ylabel('z');
    hold all;
    
    % histogram for the delay-values in the image
    ax(2) = subplot(222);
    
    hist(s2ns(aDelayData(:)),[-100:2:100]);xlim([-100 100])
    hold all;
    xlabel('ns');
    title('delay histogram');

    % mean / min / max plot
    ax(3) = subplot(223);
    if isempty(aMeanXBounds)
        xbounds = [1 size(aDelayData,2)];
    else
        xbounds = aMeanXBounds;
    end;
    
    plot(mean(s2ns(aDelayData(:,xbounds(1):xbounds(2))),2),'DisplayName','mean');
    hold all;
    plot(min(s2ns(aDelayData(:,xbounds(1):xbounds(2))),[],2),'DisplayName','min');
    plot(max(s2ns(aDelayData(:,xbounds(1):xbounds(2))),[],2),'DisplayName','max');
%     plot(median(s2ns(aDelayData(:,xbounds(1):xbounds(2))),2),'DisplayName','median');
    title('different measures across images');
%     legend(ax(3),'show','Location','SouthWest');
    if ~isempty(aDelayBounds)
        ylim(2*s2ns(aDelayBounds));
    end;
    ax(4) = subplot(224);
    
end