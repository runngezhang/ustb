function [gainA, gainB, dynA, dynB, gcA, gcB, gcPointsA, gcPointsB] = gainDynCheck(envA, envB, gainA, gainB, dynA, dynB, legA, legB, cx,cy, compressFcn,gcA,gcB, noInteraction)
    if ~exist('noInteraction','var')
        noInteraction = false;
    end;
    showProgress = true;
    if noInteraction
        showProgress = false;
    end;
    
    
    
    if exist('envB','var') && ~isempty(envB)
        assert(isequal(size(envA),size(envB)),'size of envA and envB have to be equal');
    else
        envB = envA;
    end;
    
    if ~exist('cx','var')
        cx = 5;
    end;
    if ~exist('cy','var')
        cy = 5;
    end;
    if ~exist('compressFcn','var') || isempty(compressFcn)
        compressFcn = @(x,gain,dyn)imagelog2(x,gain,dyn,256);
    end;
    
    function [nA,nB]=plotImages()
        imgA = compressFcn(envA.*repmat(10.^(gcA.'/20),1,size(envA,2)), gainA, dynA);
        imgB = compressFcn(envB.*repmat(10.^(gcB.'/20),1,size(envB,2)), gainB, dynB);
        
        nA_raw = hist(imgA(:),1:2:256+0.5);
%         nA = medfilt1(nA_raw,6);
        nB_raw = hist(imgB(:),1:2:256+0.5);
%         nB = medfilt1(nB_raw,6);
        nA = nA_raw;
        nB = nB_raw;
%         isZ = nB_raw==0;
%         isnZ = nB_raw ~= 0;
%         spikes = [0 isZ 0] & [isnZ 0 0] & [0 0 isnZ];
%         nB = nB_raw;
%         nB(spikes) = nB_raw([spikes(1:end-2)]);
%         idx = 1:256;
        if ~showProgress
            return;
        end;
        image(imgA.*single(checkBoard)+imgB.*single(~checkBoard),'Parent',plotax);
        title(plotax,sprintf('%s vs %s',legA,legB));
        if ~ishandle(Aax) || ~ishandle(Bax)
            return;
        end;
        image(imgA,'Parent',Aax);
        title(Aax,sprintf('%s',legA));
        image(imgB,'Parent',Bax);
        title(Bax,sprintf('%s',legB));
        plot(Ahax,nA);
        line(1:numel(nA),nB,'color',[1 1 1]*0.5,'parent',Ahax);
        plot(Bhax,nB);
        line(1:numel(nA),nA,'color',[1 1 1]*0.5,'parent',Bhax);
        if setylim
            ylim(Ahax,[0 max(max(nA(2:end))*1.1, max(nB(2:end))*1.1)]);
            ylim(Bhax,[0 max(max(nA(2:end))*1.1, max(nB(2:end))*1.1)]);
        end;       
    end
    function adjGain(aGainA, aGainB, aGcA, aGcB)
        gainA = aGainA;
        gainB = aGainB;
        if ~isempty(aGcA)
            gcA = aGcA;
        end;
        if ~isempty(aGcB)
            gcB = aGcB;
        end;
        plotImages();
    end
    function adjDyn(aDynA, aDynB)
        dynA = aDynA;
        dynB = aDynB;
        plotImages();
    end
    function autoAdjGain()
        [nA, nB] = plotImages();
        prct = 0.06;
        % find 15% idx
        bpidx = 128-find(cumsum(fliplr(nB))/sum(nB(:))>prct,1,'first');
        apidx = 128-find(cumsum(fliplr(nA))/sum(nA(:))>prct,1,'first');
        dg = 2;
        if bpidx > apidx
            db = -dg;
        else
            db = +dg;
        end

        improvegain = true;
        compFunc = @(xi,yi,xn,yn)abs(xi-yi);
        compFunc = @(xi,yi,xn,yn)sum(abs(xn(max(xi,yi):end)-yn(max(xi,yi):end)));
        while improvegain
            bpidx_save = bpidx;
            nB_save = nB;
            gainB = gainB + db;
            [nA, nB] = plotImages();
            bpidx = 128-find(cumsum(fliplr(nB))/sum(nB(:))>prct,1,'first');
            if (sign(db) == -1) && bpidx < apidx
                if compFunc(bpidx_save,apidx,nB_save,nA) < compFunc(bpidx, apidx, nB,nA)
                   gainB = gainB-db;
                end;
                improvegain = false;
            elseif (sign(db) == 1) && bpidx > apidx
                if compFunc(bpidx_save, apidx, nB_save,nA) < compFunc(bpidx,apidx, nB,nA)
                   gainB = gainB-db;
                end;
                improvegain = false;
            else
                pause(0.1);
            end
        end
        plotImages();
    end
    checkBFcn = @(x,y)((~mod(x.*y,2)&(mod(y,2)~=1)&(mod(x,2)==1))|((mod(x.*y,2)==0)&(mod(y,2)~=0))).';
    if ~exist('gcA','var')
      gcA = zeros(1,size(envA,1));
    end
    if ~exist('gcB','var')
      gcB = zeros(1,size(envB,1));
    end
    
    % setup checker board
    mx = size(envA,1);
    my = size(envA,2);
    [x,y] = meshgrid(1:mx,1:my);
    x = ceil(x/mx*cx);
    y = ceil(y/my*cy);
    checkBoard = checkBFcn(x,y);
    if showProgress
        figh = figure;
        scaleFigure(2,1,figh);
        pos = get(figh,'position');
        posimg = [pos(1) pos(2)-150 pos(3:4)];
        set(figh,'position',posimg);
        Aax = subplot(1,2,1);
        axis();
        Bax = subplot(1,2,2);
        axis();

        fighMeta = figure;
        scaleFigure(2,0.5,fighMeta)
        pos = get(fighMeta,'position');
        set(fighMeta,'position',[pos(1) pos(2)-530 pos(3:4)]);
        Ahax = subplot(1,7,[1:2],'Parent',fighMeta);
        Bhax = subplot(1,7,[3:4],'Parent',fighMeta);
        plotax = subplot(1,7,[5:6],'Parent',fighMeta);
        legax = subplot(1,7,7,'Parent',fighMeta);
        setylim = true;
        plotImages();
        colormap(Aax,gray(256));
        colormap(plotax,gray(256));

        title(sprintf('%s vs %s',legA,legB));

        imagesc(checkBoard,'Parent',legax);
        title(legax,sprintf('legend: - black == "%s"',legB));
    end;
    autoAdjGain();
    
    gcPointsB.x = gainLineEdit.computeXvector(8,size(gcA,1));
    gcPointsB.y = gcB(gcPointsB.x);
    gcPointsA.x = [];
    gcPointsA.y = [];
    
    if showProgress
        if ~noInteraction
            gdc = GainDynControl(@adjGain,@adjDyn,gainA, gainB, dynA, dynB,size(envA,1),gcA,gcB);
            pause(0.01)
            unitsave = get(gdc.GuiHdl,'units');
            set(gdc.GuiHdl,'units','pixels');
            pos = get(gdc.GuiHdl,'position');
            set(gdc.GuiHdl,'position',[pos(1) 900 pos(3:4)]);
            set(gdc.GuiHdl,'units',unitsave);
    
            guihdl = gdc.guiCurveR.figh;
            unitsave = get(guihdl,'units');
            set(guihdl,'units','pixels');
            pos = get(guihdl,'position');
            set(guihdl,'position',[posimg(1)+posimg(3) posimg(2) pos(3:4)]);
            set(guihdl,'units',unitsave);
    
            gdc.setGainLock(true);
            gdc.setDynLock(true);
            waitfor(gdc.GuiHdl);
            gcPointsB.x = gdc.guiCurveR.cPtsX;
            gcPointsB.y = gdc.guiCurveR.cPtsY;
            gcPointsA.x = [];
            gcPointsA.y = [];
            close(gdc.guiCurveR.figh);
        end;
        pause(0.2);
        close(figh);
        close(fighMeta);
    end;
end