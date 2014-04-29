function status = fig2tikz(varargin)

% status = fig2tikz(hfig,filename)
%   Print figure with handle 'hfig' to pgf-file 'filename'
% status = fig2tikz(filename)
%   Print current figure to pgf-file 'filename'
% status = fig2tikz(hfig)
%   Print figure to screen
% status = fig2tikz
%   Print current figure to screen
% status = fig2tikz(hax,filename)
%   Print current axes to pgf-file 'filename'
%
% File extension '.pgf' will be added if no extension is given.
%
% return argument
%   status -- 0 indicates failure, 1 success
%
% Known issues:
% - Cannot handle 3D plots (look into Asymptote to do this)
%
% Uses plotboxpos by Kelly Kearney
%

% Last updated 20011-09-19, Øyvind Standal

status = 0;

%% Variables common for all nested functions (like global)
error(nargchk(0,2,nargin,'struct'));
if nargin && any(ishandle(varargin{1})),
    [h,varargin] = deal(varargin{1},varargin(2:end));
else
    h = gcf;
end
if ~isempty(varargin),
    filename = varargin{1};
    [fpath,fname,ext] = fileparts(filename);
    if isempty(ext),
        filename = fullfile(fpath,[fname,'.pgf']);
    end
    fid = fopen(filename,'w+');
else
    fid = 1; % std output (screen)
end

if strcmpi(get(h,'type'),'axes'),
    h_fig = ancestor(h,'figure');
    h_child = h; % draw only this one child axis
elseif strcmpi(get(h,'type'),'figure'),
    h_fig = h;%copyobj(h,0);
    %set(h_fig,'Visible','off');
    h_child = allchild(h_fig); % draw all child axes
end

% suspend and store properties of figure
figstate0 = uisuspend(h_fig);

Colormap = get(h_fig,'Colormap');
scalefunc = containers.Map({'linear','log'},{@(x)x , @log10});
rescale = @(x,y)feval(scalefunc(x),y);
% Declare variables that are to be available for all nested functions
[P0,Width,Height,XLim,YLim,CLim,XScale,YScale,XScaleFactor,YScaleFactor,XDir,YDir] = deal([]);
colordict = {'white',[1,1,1]; ...
             'black',[0,0,0]; ...
             'red'  ,[1,0,0]; ...
             'green',[0,1,0]; ...
             'blue' ,[0,0,1]};
maxploterror = 0.2; % in fraction of linewidth
%%

% Define layers that allow us to put certain drawings on top of others,
% regardless of the order in which they are drawn.
fprintf(fid,'%% Created %s by fig2tikz.m\n',datestr(now,31));
fprintf(fid,'%% Øyvind Standal\n\n');
fprintf(fid,'%s\n','\pgfdeclarelayer{bottom}');
fprintf(fid,'%s\n','\pgfdeclarelayer{top}');
fprintf(fid,'%s\n\n','\pgfsetlayers{bottom,main,top}');

fprintf(fid,'%s\n\n','\usepgflibrary{plotmarks}'); % load plot marks

fprintf(fid,'%s\n','\def\tikzscale{0.95}');
fprintf(fid,'%s\n\n','\begin{tikzpicture}[scale=\tikzscale] % change scale of picture if desired');

try
    draw_handles(h_child);
catch err,
    if isempty(intersect(fid,0:2)),
        fclose(fid);
    end
    err.rethrow();
end

fprintf(fid,'%s\n\n','\end{tikzpicture}');
if isempty(intersect(fid,0:2)),
    fclose(fid);
end

% restore figure state
uirestore(figstate0);

status = 1;

%%
    function draw_handles(handles)
        handles = handle(handles);
        % reverse order, latest objects should be on top (except as defined
        % by \pgfonlayer)
        for j = numel(handles):-1:1,
            if ishandle(handles(j)),
                switch handles(j).Type,
                    case 'axes',
                        draw_axes(handles(j));
                    case 'line',
                        draw_line(handles(j));
                    case 'image',
                        draw_image(handles(j));
                    case 'patch',
                        draw_patch(handles(j));
                    case 'text',
                        draw_text(handles(j));
                    otherwise
                        % Do nothing
                end
                draw_handles(allchild(handles(j)));
            end
        end
    end % draw_handles

%%
    function draw_axes(hax)

        pos = plotboxpos(hax,'centimeters');
        P0 = pos(1:2);
        Width = pos(3);
        Height = pos(4);
        avgridstyles = {'-','solid';'--','dashed';':','dotted'}; % more are possible
        gridstyle = avgridstyles{strcmpi(hax.Grid,avgridstyles(:,1)),2};
        minorgridstyle = avgridstyles{strcmpi(hax.MinorGridLineStyle,avgridstyles(:,1)),2};
        
        lw = hax.LineWidth;
        tl = 8 * lw; % tick length
        mtl = tl/4; % minor tick length

        XScale = hax.XScale;
        YScale = hax.YScale;
        XLim = rescale(XScale,hax.XLim);
        YLim = rescale(YScale,hax.YLim);
        CLim = get(hax,'CLim');

        XDir = hax.XDir;
        YDir = hax.YDir;
        
        % important: colors must always be defined in main pgflayer
        axColor = getColorName(get(hax,'Color'));
        XColor = getColorName(get(hax,'XColor'));
        YColor = getColorName(get(hax,'YColor'));

        XScaleFactor = diff(XLim)/Width;
        xtick = (rescale(XScale,hax.XTick)-XLim(1)) / XScaleFactor;
        xticklabel = strtrim(cellstr(hax.XTickLabel));
        if strcmpi(XScale,'log'),
            xticklabel = cellfun(@(x)sprintf('$10^{%s}$',x),xticklabel,'un',0);
        end
        xticklabel = cellfun(@txt2latex,xticklabel,'uniformoutput',false);
        if strcmpi(XDir,'reverse'),
            xtick = Width - xtick;
        end
        try % should find better solution
            X = merge(xtick,xticklabel);
            X = deblank(sprintf('%.3f/%s, ',X{:}));
            X = X(1:end-1); % drop last comma
        catch
            X = '';
        end

        YScaleFactor = diff(YLim)/Height;
        ytick = (rescale(YScale,hax.YTick)-YLim(1)) / YScaleFactor;
        yticklabel = strtrim(cellstr(hax.YTickLabel));
        if strcmpi(YScale,'log'),
            yticklabel = cellfun(@(x)sprintf('$10^{%s}$',x),yticklabel,'un',0);
        end
        yticklabel = cellfun(@txt2latex,yticklabel,'uniformoutput',false);
        if strcmpi(YDir,'reverse'),
            ytick = Height - ytick;
        end
        try
            Y = merge(ytick,yticklabel);
            Y = deblank(sprintf('%.3f/%s, ',Y{:}));
            Y = Y(1:end-1); % drop last comma
        catch
            Y = '';
        end

        if isa(handle(hax),'scribe.legend'),
            fprintf(fid,'%s\n','% Draw legend');
            axdrawcmd = '\draw[fill=white,draw=black,line width=\lw] (P0) rectangle +(1.1*\width,\height);';
            X = []; Y = []; % sometimes there are dummy ticks in legend
            layer = 'top'; % place on top
        else
            fprintf(fid,'%s\n','% Draw axes');
            XAxColName = getColorName(hax.XColor);
            YAxColName = getColorName(hax.YColor);
            if strcmpi(hax.Box,'on'),
                axdrawcmd = sprintf(['\\draw[draw=%s,line width=\\lw] [shift={(P0)}] (0,0) -- +(\\width,0) (0,\\height) -- +(\\width,0);\n', ...
                                         '\\draw[draw=%s,line width=\\lw] [shift={(P0)}] (0,0) -- +(0,\\height) (\\width,0) -- +(0,\\height);'], ...
                                         XAxColName,YAxColName);
            else
                if strcmpi(hax.XAxisLocation,'bottom'),
                    xaxdrawcmd = sprintf('\\draw[draw=%s,line width=\\lw] (P0) -- +(\\width,0);',XAxColName);
                else
                    xaxdrawcmd = sprintf('\\draw[draw=%s,line width=\\lw] (P0) ++(0,\\height) -- +(\\width,0);',XAxColName);
                end
                if strcmpi(hax.YAxisLocation,'left'),
                    yaxdrawcmd = sprintf('\\draw[draw=%s,line width=\\lw] (P0) -- +(0,\\height);',YAxColName);
                else
                    yaxdrawcmd = sprintf('\\draw[draw=%s,line width=\\lw] (P0) ++(\\width,0) -- +(0,\\height);',YAxColName);
                end
                axdrawcmd = sprintf('%s\n%s',xaxdrawcmd,yaxdrawcmd);
            end
            layer = lower(get(hax,'Layer'));
        end
        %fprintf(fid,'%s\n\n',\begin{scope}');
        fprintf(fid,'%s\n','% Define some constants');
        fprintf(fid,'\\def \\width {%.3f};\n',Width);
        fprintf(fid,'\\def \\height {%.3f};\n',Height);
        fprintf(fid,'\\def \\lw {%gpt}; %% line width\n',lw);
        fprintf(fid,'\\def \\tl {%gpt}; %% tick length\n',tl);
        fprintf(fid,'\\def \\mtl {%gpt}; %% minor tick length\n',mtl);
        % fprintf(fid,'\\definecolor{%s}{rgb}{%.3f,%.3f,%.3f};\n',getColorName(),color);
        fprintf(fid,'\\coordinate (P0) at (%.3f,%.3f); %% define axis origin\n\n',P0);
        
        if strcmpi(hax.Visible,'off'),
            % axis is not drawn, but keep defined constants
            return
        end
        
        function mticks = getminorticks(ticks,lims,scale)
            if numel(ticks) < 2,
                mticks = [];
                return
            end
            subgrid = [2,5,10];
            if nargin < 2 || strcmpi(scale,'linear'),
                dt = diff(ticks(1:2));
                nsubs = subgrid(find(dt/0.25>=subgrid ,1,'last'));
                if isempty(nsubs), mticks = []; return; end
                dd = dt/nsubs;
                mticks = unique([ticks(1):-dd:lims(1),ticks(1):dd:lims(2)]);
                mticks = setdiff(mticks,lims);
            else % log
                dt = diff(ticks);
                ticks = [ticks(1)-dt(1),ticks,ticks(end)+dt(end)];
                mticks = unique(bsxfun(@plus,ticks(1:end-1)',diff(ticks)'*log10(1:10)));
                mticks = mticks(mticks>lims(1) & mticks<lims(end));
            end
        end
            
        % fprintf(fid,'%s\n','\begin{pgfonlayer}{bottom}');
        % fprintf(fid,'%s\n',axdrawcmd);
        % fprintf(fid,'%s\n','\end{pgfonlayer};');
        
        fprintf(fid,'\\begin{pgfonlayer}{%s}\n',layer);
        if ~isempty(X) || ~isempty(Y),
            fprintf(fid,'%s\n','% Draw tick marks');
        end
        % minor grid/tick
        xminortick = getminorticks(xtick,[0,diff(XLim)]/XScaleFactor,hax.XScale);
        if ~isempty(xminortick),
            doDrawMinor = [0,0];
            minorstr = sprintf('%.3f,',xminortick);
            minorstr = sprintf('\\foreach \\x in {%s}{',minorstr(1:end-1));
            if strcmpi(hax.XMinorTick,'on'),
                doDrawMinor(1) = 1;
                minorstr = sprintf('%s\n\t\\draw[%s,line width=\\lw] (P0) ++(\\x,0) -- +(0,\\mtl);',minorstr,XColor);
            end
            if strcmpi(hax.XMinorGrid,'on'),
                doDrawMinor(2) = 1;
                minorstr = sprintf('%s\n\t\\draw[gray,thin,%s] (P0) ++(\\x,0) -- +(0,\\height);',minorstr,minorgridstyle);
            end
            if any(doDrawMinor),
                fprintf(fid,'%s\n\t}\n\n',minorstr);
            end
        end
        % minor grid/tick
        yminortick = getminorticks(ytick,[0,diff(YLim)]/YScaleFactor,hax.YScale);
        if ~isempty(yminortick),
            doDrawMinor = [0,0];
            minorstr = sprintf('%.3f,',yminortick);
            minorstr = sprintf('\\foreach \\y in {%s}{',minorstr(1:end-1));
            if strcmpi(hax.YMinorTick,'on'),
                doDrawMinor(1) = 1;
                minorstr = sprintf('%s\n\t\\draw[%s,line width=\\lw] (P0) ++(0, \\y) -- +(\\mtl, 0);',minorstr,YColor);
            end
            if strcmpi(hax.YMinorGrid,'on'),
                doDrawMinor(2) = 1;
                minorstr = sprintf('%s\n\t\\draw[gray,thin,%s] (P0) ++(0,\\y) -- +(\\width,0);',minorstr,minorgridstyle);
            end
            if any(doDrawMinor),
                fprintf(fid,'%s\n\t}\n\n',minorstr);
            end
        end
        if ~isempty(X),
            fprintf(fid,'\\foreach \\x/\\xlabel in {%s}{\n',X);
            if strcmpi(hax.XAxisLocation,'bottom'),
                fprintf(fid,'\t\\draw[%s,line width=\\lw] (P0) ++(\\x, .5*\\tl) -- +(0, -\\tl) node[anchor=north] {\\xlabel};',XColor);
            else
                fprintf(fid,'\t\\draw[%s,line width=\\lw] (P0) ++(0,\\height) +(\\x,-.5*\\tl) -- +(\\x, .5*\\tl) node[anchor=south] {\\xlabel};',XColor);
            end
            if any(strcmpi({hax.XGrid,hax.XMinorGrid},'on')),
                fprintf(fid,'\n\t\\draw[black,thin,%s] (P0) ++(\\x,0) -- +(0,\\height);\n\t}\n',gridstyle);
            else
                fprintf(fid,'\n\t}\n');
            end
        end
        
        if ~isempty(Y),
            fprintf(fid,'\\foreach \\y/\\ylabel in {%s}{\n',Y);
            if strcmpi(get(hax,'YAxisLocation'),'left'),
                fprintf(fid,'\t\\draw[%s,line width=\\lw] (P0) ++(.5*\\tl, \\y) -- +(-\\tl, 0) node[anchor=east] {\\ylabel};',YColor);
            else
                fprintf(fid,'\t\\draw[%s,line width=\\lw] (P0) ++(\\width,0) +(-.5*\\tl, \\y) -- +(0.5*\\tl, \\y) node[anchor=west] {\\ylabel};',YColor);
            end
            if any(strcmpi({hax.YGrid,hax.YMinorGrid},'on')),
                fprintf(fid,'\n\t\\draw[black,thin,%s] (P0) ++(0,\\y) -- +(\\width,0);\n\t}\n',gridstyle);
            else
                fprintf(fid,'\n\t}\n\n');
            end 
        end
        fprintf(fid,'%s\n',axdrawcmd);
        fprintf(fid,'%s\n\n','\end{pgfonlayer}');
    end % draw_axes


%% Image
    function draw_image(hi)

        persistent fileno
        
        if isempty(fileno),
            fileno = 1;
        end
        
        cdata = get(hi,'CData');
        map = get(hi,'CDataMapping');
        
        xdata = get(hi,'xdata');
        if size(cdata,2)==1,
            dx = XScaleFactor * Width;
            xdata = mean(xdata);
            xin = 1;
        else
            dx = diff(xdata([1,end]))/(size(cdata,2)-1);
            xdata = linspace(xdata(1),xdata(end),size(cdata,2));
            xin = xdata+dx/2 > XLim(1) & xdata-dx/2 < XLim(2);
        end
        ydata = get(hi,'ydata');
        if size(cdata,1)==1,
            dy = Height/2;
            yin = 1;
        else
            dy = diff(ydata([1,end]))/(size(cdata,1)-1);
            ydata = linspace(ydata(1),ydata(end),size(cdata,1));
            yin = ydata+dy/2 > YLim(1) & ydata-dy/2 < YLim(2);
        end
        
        x = xdata(xin);
        y = ydata(yin);
        c = cdata(yin,xin);
        
        if size(c,3)==3,
            rgb = c; % haven't checked if this works
        elseif strcmpi(map,'direct'),
            rgb = Colormap(c(:),:);
        else % strcmpi(map,'scaled'),
            rgb = Colormap(round( interp1(CLim,[1,size(Colormap,1)],c(:)) ),:);
        end
        
        if strcmpi(XDir,'reverse'),
            x = -x + XLim(2);
        else
            x = x - XLim(1);
        end
        if strcmpi(YDir,'reverse'),
            y = -y + YLim(2);
        else
            y = y - YLim(1);
        end
        
        if numel(cdata) < 0, % not too big image
            [X,Y] = meshgrid(x,y);
            V = [X(:)'/XScaleFactor;Y(:)'/YScaleFactor;rgb'];
            imgstr = deblank(sprintf('%.3f/%.3f/{%.3f,%.3f,%.3f}, ',V(:)'));
            
            fprintf(fid,'%s\n','% Draw image');
            fprintf(fid,'%s\n','\begin{scope}');
            fprintf(fid,'%s\n','\clip (P0) rectangle +(\width,\height);');
            % fprintf(fid,'\\path(%.3f,%.3f) coordinate (P0);\n',P0);
            fprintf(fid,'\\def \\dx {%.3f};\n',dx/XScaleFactor);
            fprintf(fid,'\\def \\dy {%.3f};\n',dy/YScaleFactor);
            fprintf(fid,'%s\n','\def \rectpath {+(-.5*\dx,-.5*\dy) rectangle +(.5*\dx,.5*\dy)};');
            
            fprintf(fid,'\\foreach \\x/\\y/\\rgb in {%s}\n',imgstr(1:end-1));% remove last comma
            fprintf(fid,'%s\n','{'); % loop delimiter
            % don't use getColorName here
            fprintf(fid,'\t%s\n','\definecolor{imgCol}{rgb}{\rgb};');
            fprintf(fid,'\t%s\n','\fill[color=imgCol] (P0) ++(\x,\y) \rectpath;');
            fprintf(fid,'%s\n','};'); % end loop
            fprintf(fid,'%s\n\n','\end{scope}');
        
        else % large image
            pngfilename = regexprep(filename,'\.\w+$',sprintf('%d.png',fileno));
            while exist(pngfilename,'file'),
                fileno = fileno + 1;
                pngfilename = regexprep(filename,'\.\w+$',sprintf('%d.png',fileno));
            end
            [~,pngfile_includename,pext] = fileparts(pngfilename); % ditch path before include
            rgb = reshape(rgb,numel(y),numel(x),3);
            if strcmpi(YDir,'normal'),
                rgb = flipdim(rgb,1);
            end
            if strcmpi(XDir,'reverse'),
                rgb = flipdim(rgb,2); % not sure!
            end
            imwrite(rgb,pngfilename,'png');
            fprintf(fid,'%s\n','% Draw image');
            fprintf(fid,'%s\n','\begin{scope}');
            fprintf(fid,'%s\n','\clip (P0) rectangle +(\width,\height);');
            fprintf(fid,'%s\n','\pgfmathsetmacro{\pngwidth}{\tikzscale*\width}');
            fprintf(fid,'%s\n','\pgfmathsetmacro{\pngheight}{\tikzscale*\height}');
            fprintf(fid,'\\node [anchor=south west, inner sep=0pt] at (P0) {\\includegraphics[width=\\pngwidth cm,height=\\pngheight cm]{%s}};\n',[pngfile_includename,pext]);
            fprintf(fid,'%s\n\n','\end{scope}');
        end
    end

%% Patch
    function draw_patch(h,xdata,ydata,cdata)

        if nargin < 4,
            cdata = h.CData;
            % hack because matlab sometimes seem to give wrong shape to cdata
            cdata = shiftdim(squeeze(cdata),ndims(squeeze(cdata))-ndims(cdata));
            if ndims(cdata) < 3, % Change indexed color to truecolor rgb
                cdata = ind2rgb(cdata,Colormap,CLim);
            end
            xdata = h.XData;
            ydata = h.YData;
        end
        if size(xdata,2) > 1,
            for k = 1:size(xdata,2),
                if ~isempty(cdata),
                    draw_patch(h,xdata(:,k),ydata(:,k),cdata(:,min(end,k),:));
                else
                    draw_patch(h,xdata(:,k),ydata(:,k),[]);
                end
            end
            return % important to remember this...
        end
        reject = (isnan(xdata) | isnan(ydata));
        xdata = xdata(~reject);
        ydata = ydata(~reject);
        %cdata = cdata(~reject);
        if max(xdata) < XLim(1) || min(xdata) > XLim(2) || ...
                max(ydata) < YLim(1) || min(ydata) > YLim(2),
            % nothing to draw here
            return
        end
        % check if clipping is necessary
        doClip = min(xdata) < XLim(1) || max(xdata) > XLim(2) || ...
                min(ydata) < YLim(1) || max(ydata) > YLim(2);
        if doClip,
            [xdata,ydata] = clip_patch(xdata,XLim,ydata,YLim);
        end
            
        edgecolor = h.EdgeColor;
        facecolor = h.FaceColor;

        x = (xdata - XLim(1)) / XScaleFactor;
        y = (ydata - YLim(1)) / YScaleFactor;
        
        [x,y] = prune_line(x,y,maxploterror*h.LineWidth*2.54/72);
        if x(1)==x(end) && y(1)==y(end), % if path is closed, open it again
            x = x(1:end-1); y = y(1:end-1);
        end
        if numel(x) < 3 || numel(y) < 3,
            return
        end
        if strcmpi(XDir,'reverse'),
            x = diff(XLim)/XScaleFactor - x;
        end
        if strcmpi(YDir,'reverse'),
            y = diff(YLim)/YScaleFactor - y;
        end

        patchpath = [sprintf('(%.3f,%.3f) -- ',[x(:)';y(:)']),'cycle'];

        if strcmp(edgecolor, 'flat'),
            edgecolor = squeeze(cdata(1,1,:));
        end

        if strcmp(facecolor, 'flat'),
            facecolor = squeeze(cdata(1,1,:));
        end

        fill = ''; draw = '';
        if isnumeric(edgecolor),
            % fprintf(fid,'\\definecolor{edgeCol}{rgb}{%.3f,%.3f,%.3f};\n',edgecolor);
            edgecolname = getColorName(edgecolor);
            draw = sprintf('[draw=%s,line width=%.3f]',edgecolname,h.LineWidth);
        end
        if isnumeric(facecolor),
            % fprintf(fid,'\\definecolor{faceCol}{rgb}{%.3f,%.3f,%.3f};\n',facecolor);
            facecolname = getColorName(facecolor);
            fill = sprintf('[fill=%s]',facecolname);
        end
        fprintf(fid,'%s\n','% Draw patch');
        if 0,%doClip,
            fprintf(fid,'%s\n','\begin{scope}');
            fprintf(fid,'%s\n','\clip (P0) rectangle +(\width,\height);');
        end
        fprintf(fid,'\\path[shift={(P0)}] %s %s {%s};\n',draw,fill,patchpath);
        if 0,%doClip,
            fprintf(fid,'%s\n','\end{scope}');
        end
        fprintf(fid,'\n');

    end

%% Text
    function draw_text(ht)

        str = cellstr(deblank(get(ht,'String')));
        doSplit = ( numel(str) > 1 ); % Split string across several lines
        if isempty(char(str)) || strcmpi(get(ht,'Visible'),'off'),
            return
        end
        str = deblank(sprintf('%s ',str{:}));
        str = txt2latex(str);
        
        pos = getprop(ht,'position','centimeters');
        extent = 1.5 * getprop(ht,'Extent','Centimeters');
        rotation = get(ht,'Rotation');
        valign = get(ht,'VerticalAlignment');
        halign = get(ht,'HorizontalAlignment');

        vdict = {'bottom','south';'baseline','south';'middle','';'cap','north';'top','north'};
        hdict = {'left','west';'center','';'right','east'};
        anchor = strtrim(sprintf('%s %s',vdict{strcmpi(valign,vdict),2}, hdict{strcmpi(halign,hdict),2}));
        if isempty(anchor),
            anchor = 'base'; % default
        end

        % Draw
        fprintf(fid,'%s\n','% Draw text');
        setLayer = isa(handle(get(ht,'Parent')),'scribe.legend');
        if setLayer,
            fprintf(fid,'%s\n','\begin{pgfonlayer}{top}');
        end
        txtcolor = get(ht,'Color');
        txtcolname = getColorName(txtcolor);
        fprintf(fid,'\\path (P0) +(%.3f,%.3f) coordinate (X0); %% define axis origin\n',pos(1:2));
        % fprintf(fid,'\\definecolor{myCol}{rgb}{%.3f,%.3f,%.3f};\n',color);
        if doSplit, % possibly split string into several lines
            fprintf(fid,'\\node[%s,anchor=%s,rotate=%.3f,text width=%.0fcm,minimum height=%.0f] at (X0) {%s};\n', ...
                txtcolname,anchor,rotation,extent(3),extent(4),str);
        else % everything on one line
            fprintf(fid,'\\node[%s,anchor=%s,rotate=%.3f] at (X0) {%s};\n',txtcolname,anchor,rotation,str);
        end
        if setLayer,
            fprintf(fid,'%s\n','\end{pgfonlayer}');
        end
        fprintf(fid,'\n');
    end


%% Line plotting
    function draw_line(h,xdata,ydata)

        if ~strcmpi(h.Visible,'on') || ( strcmpi(h.LineStyle,'none') && strcmpi(h.Marker,'none') ),
            return
        end
        linewidth = h.LineWidth/1.0;
        edgeCol = h.MarkerEdgeColor;
        faceCol = h.MarkerFaceColor;
        if strcmpi(edgeCol,'none'),
            edgeCol = '';
        elseif strcmpi(edgeCol,'auto'),
            edgeCol = h.Color;
        end
        if strcmpi(faceCol,'none'),
            faceCol = '';
        elseif strcmpi(faceCol,'auto'),
            faceCol = h.Color;
        end
        
        if nargin < 3,
            xdata = rescale(XScale,h.XData);
            ydata = rescale(YScale,h.YData);
        end
        [x,y,inside_index] = clip(xdata,XLim,ydata,YLim);
        if iscell(x),
            for j = 1:length(x),
                % if inside_index{j} etlerannet
                draw_line(h,x{j},y{j});
            end
            return
        end
        
        if isempty(x) || isempty(y),
            % nothing to draw here
            return
        end
        doClip = ~all(inside_index); % if any datapoints are outside axes, clip

        x = (x - XLim(1)) / XScaleFactor;
        y = (y - YLim(1)) / YScaleFactor;
        
        if strcmpi(h.Marker,'none'),
            % Cannot remove data when markers are present
            [x,y] = prune_line(x,y,maxploterror*linewidth*2.54/72);
        end
        data = [x(:)';y(:)'];
        datastr = deblank(sprintf('(%.3f,%.3f) ',data(:)));

        markerstr = getMarkerString(h.Marker,h.MarkerSize,getColorName(edgeCol),getColorName(faceCol));
        linestylestr = getLineStyleString(h.LineStyle);

        fprintf(fid,'%s\n','% Draw line');
        setLayer = isa(handle(ancestor(h,'axes')),'scribe.legend');
        if setLayer,
            fprintf(fid,'%s\n','\begin{pgfonlayer}{top}');
        end
        lineColName = getColorName(h.Color);
        if doClip,
            fprintf(fid,'%s\n','\begin{scope}');
            fprintf(fid,'%s\n','\clip (P0) rectangle +(\width,\height);');
        end
        if isempty(linestylestr), % only marks, don't draw line
            fprintf(fid,'\\def \\myPlot {plot[only marks] coordinates {%s}};\n',datastr);
        else
            fprintf(fid,'\\def \\myPlot {plot coordinates {%s}};\n',datastr);
        end
        fprintf(fid,'\\draw[%s,color=%s,shift={(P0)},line width=%.3f,mark=%s] \\myPlot;\n', ...
            linestylestr,lineColName,linewidth,markerstr);
        if doClip,
            fprintf(fid,'%s\n','\end{scope}');
        end
        if setLayer,
            fprintf(fid,'%s\n','\end{pgfonlayer}');
        end
        fprintf(fid,'\n');   
        
    end

%% Color management
    function colname = getColorName(rgb,addtodict)
        % If color rgb has already been used, find its given name from colordict.
        % If not, add color to colordict, and generate name for color.
        % If addtodict is false, do not add color to colordict (e.g. when color is defined inside a scope)
        addtodict = nargin < 2 || addtodict;
        if isempty(rgb) || any(isnan(rgb)) || strcmpi(rgb,'none'),
            colname = 'none';
            return
        end
        rgb = rgb(:)';
        ind = find(cellfun(@(x)all(x==rgb),colordict(:,2)));
        if ~isempty(ind),
            colname = colordict{ind,1};
        else
            colname = sprintf('myCol%d',max(1,size(colordict,1)-4));
            if addtodict,
                colordict(end+1,:) = {colname,rgb};
            end
            % define color in pgf-file
            fprintf(fid,'\\definecolor{%s}{rgb}{%.3f,%.3f,%.3f};\n',colname,rgb);
        end
    end


end % main


%% Subfunctions
function texstr = txt2latex(str)
% Add $-signs where appropriate to get latex-string
% Only deals with numbers, super- and subscript

% Don't try to mess with a string that's already in latex format or containing '\_' stuff
ignore_strings = {'\$','\_'};
if ~all(cellfun(@isempty,regexp(str,ignore_strings, 'once' ))),
    texstr = str;
    return
end
% texstr = regexprep(str,'(\\\S+)|(\S+[_^](([^{])|({.+})))|\<(\d+(\.\d+)?)\>','$$1$');
str = regexprep(str,'(?<!\\)\%','\\%'); % replace '%' by '\%'
texstr = regexprep(str,'(\S+[_^](([^{])|({.+})))|\<(-?\d+(\.\d+)?)\>','$$1$');
end

function rgb = ind2rgb(x,cmap,clim)
rgb = NaN(length(x),3);
nanind = isnan(x);
x = x(~nanind);
if nargin==3, % scaled indexing
    x = reshape(round( interp1(clim,[1,length(cmap)],x(:)) ),size(x));
end
rgb(~nanind,:) = reshape(cmap(x(:),:),[size(x),3]);
end

function str = sprintf(varargin)
str = builtin('sprintf',varargin{:});
str = regexprep(regexprep(str,'(\d+\.\d+?)0*\>','$1'),'\.0\>','');
end

function fprintf(fid,varargin)
str = sprintf(varargin{:});
builtin('fprintf',fid,'%s',str);
end

function C = merge(varargin)

if nargin==1,
    C = varargin{1};
    return
end

if all(cellfun(@isnumeric,varargin)),
    C = zeros(1,nargin*numel(varargin{1}));
    for j = 1:nargin,
        C(j:nargin:end) = varargin{j};
    end
    return
elseif all(cellfun(@ischar,varargin)),
    C = char(varargin{:});
    C = C(:)';
    return
end

C = cell(1,nargin*length(varargin{1}));
v = cell(1,nargin);
for j = 1:nargin,
    if isempty(varargin{j}),
        continue
    end
    if isnumeric(varargin{j}),
        v{j} = num2cell(varargin{j});
    elseif ischar(varargin{j}),
        v{j} = cellstr(varargin{j});
    else
        v{j} = varargin{j};
    end
    C(j:nargin:end) = v{j};
end

if all(cellfun(@isnumeric,varargin)),
    C = cell2mat(C);
elseif all(cellfun(@ischar,varargin)),
    C = char(C);
end
end


function [x,y] = prune_line(x,y,tol)
% Remove datapoints that can be replaced by a straight line.

if nargin < 3,
    tol = 0.001;
end
if numel(x) < 3 || numel(y) < 3,
    return
end
x = x(:); y = y(:);
ind = 0;
while ~all(ind),
    ind = (ldist(x,y)>tol);
    ind(2:2:end) = true; % remove at most every other sample
    x = x([true;ind(:);true]);
    y = y([true;ind(:);true]);
end
% Final pass over coordinates
ind = (ldist(x,y)>tol);
ind(1:2:end) = true; % remove at most every other sample
x = x([true;ind(:);true]);
y = y([true;ind(:);true]);
    function val = ldist(x,y)
        val = abs((x(3:end)-x(1:end-2)).*(y(1:end-2)-y(2:end-1)) - ...
            (x(1:end-2)-x(2:end-1)).*(y(3:end)-y(1:end-2))) ./ ...
            sqrt((x(3:end)-x(1:end-2)).^2 + (y(3:end)-y(1:end-2)).^2);
        val(isnan(val)) = inf;
    end
end


function [xout,yout] = clip_patch(x,xlim,y,ylim)
xout = min(max(x,xlim(1)),xlim(2));
yout = min(max(y,ylim(1)),ylim(2));
end

%%
function [xout,yout,in] = clip(x,xlim,y,ylim)
% Trim data that is outside xlim and ylim, also NaN and Inf
x = x(:)'; y = y(:)';
if isempty(x) || isempty(y),
    [xout,yout,in] = deal([]);
    return
end
isvalid = ~isnan(x) & ~isinf(x) & ~isnan(y) & ~isinf(y);
s1 = find([isvalid(1),diff(isvalid)]==1);
s2 = find([~isvalid(1),diff(~isvalid)]==1);
if isempty(s1),
    [xout,yout,in] = deal([]);
    return
elseif length(s1) > 1,
    if s2(1) < s1(1),
        s2 = s2(2:end);
    end
    if s1(end) > s2(end),
        s2(end+1) = length(isvalid)+1;
    end
    [xout,yout,in] = deal(cell(size(s1)));
    for j = 1:length(s1),
        ind = s1(j):s2(j)-1;
        [xout{j},yout{j},in{j}] = clip(x(ind),xlim,y(ind),ylim);
    end
else
    x = x(isvalid);
    y = y(isvalid);
    % Remove points outside axis
    in = inpolygon(x,y,xlim([1,2,2,1,1]),ylim([1,1,2,2,1]));
    keep = true(size(x));
    for j = 2:length(in)-1,
        xj = x(j:j+1); yj = y(j:j+1);
        keep(j) = any(in(j-1:j+1)) || ...
            line_cross(xj,yj,xlim([1,2]),ylim([1,1])) || ...
            line_cross(xj,yj,xlim([2,2]),ylim([1,2])) || ...
            line_cross(xj,yj,xlim([2,1]),ylim([2,2])) || ...
            line_cross(xj,yj,xlim([1,1]),ylim([2,1]));
    end
    if length(keep) > 1,
        keep([1,end]) = keep([2,end-1]);
    end
    xout = x(keep);
    yout = y(keep);
end

    function out = line_cross(a,b,c,d)
        % Logical. True if line a-b crosses line c-d.
        out = false;
        M = [-diff(a),diff(c);-diff(b),diff(d)];
        if rank(M) == 2,
            v = [a(1)-c(1);b(1)-d(1)];
            t = M\v;
            out = all(t >=0 & t <= 1);
        end
    end

end


%% Subfunctions
function str = getLineStyleString(linestyle)
switch lower(linestyle),
    case 'none',
        str = '';
    case '-',
        str = 'solid';
    case ':',
        str = 'dotted';
    case '--',
        str = 'dashed';
    case '-.',
        str = 'dash pattern=on 0.8pt off 2pt on 3pt off 3pt';
    otherwise
        warning('fig2tikz:getMarkerString:unsupportedLinestyle','Unknown linestyle: %s. Using solid line.',linestyle);
end
end


function str = getMarkerString(marker,markersize,color,facecolor)
markersize = markersize/2;
switch lower(marker),
    case 'none',
        str = '';
    case '.',
        str = sprintf('*,mark size=%gpt',markersize/3);
    case '*',
        str = sprintf('asterisk,mark size=%gpt',markersize);
    case '^',
        str = sprintf('triangle,mark size=%gpt,mark options={solid,draw=%s,fill=%s}',markersize,color,facecolor);
    case '<',
        str = sprintf('triangle,mark size=%gpt,mark options={solid,draw=%s,fill=%s,rotate=90}',markersize,color,facecolor);
    case 'v',
        str = sprintf('triangle,mark size=%gpt,mark options={solid,draw=%s,fill=%s,rotate=180}',markersize,color,facecolor);
    case '>',
        str = sprintf('triangle,mark size=%gpt,mark options={solid,draw=%s,fill=%s,rotate=270}',markersize,color,facecolor);
    case {'s','square'},
        str = sprintf('square,mark size=%gpt,mark options={solid,draw=%s,fill=%s}',markersize/1.4,color,facecolor);
    case {'d','diamond'},
        str = sprintf('diamond,mark size=%gpt,mark options={solid,draw=%s,fill=%s}',markersize,color,facecolor);
    case {'p','pentagram'},
        str = sprintf('star,mark size=%gpt,mark options={solid,draw=%s,fill=%s}',markersize,color,facecolor);
    case '+',
        str = sprintf('+,mark size=%gpt,mark options={solid,draw=%s,fill=%s}',markersize,color,facecolor);
    case 'x',
        str = sprintf('x,mark size=%gpt,mark options={solid,draw=%s,fill=%s}',markersize,color,facecolor);
    case 'o',
        str = sprintf('o,mark size=%gpt,mark options={solid,draw=%s,fill=%s}',markersize,color,facecolor);
    otherwise
        warning('fig2tikz:getMarkerString:unsupportedMarker','Unknown marker: %s. Not using any marker.',marker);
        str = '';
end
if ~isempty(facecolor),
    str = regexprep(str,'^(\w+),','$1*,','once');
end
end

function pos = plotboxpos(h,units)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2006 Kelly Kearney

% Made a few minor changes (2008-08-04 OKVS)
% - Now returns position in desired units.
% - No longer leaves a dummy figure onscreen

if nargin < 2,
    units = get(h,'units');
end

axisPos = getprop(h,'Position','Pixels');
darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    dx = diff(get(h, 'XLim'));
    dy = diff(get(h, 'YLim'));
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

ftemp = figure('visible','off');
try
    temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off');
    pos = getprop(temp,'Position',units);
catch
    pos = [];
end
delete(ftemp);

end


%%
function p = getprop(h,prop,units)
units0 = get(h,'units');
set(h,'units',units);
p = get(h,prop);
set(h,'units',units0);
end
