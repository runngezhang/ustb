function make_movie2_seq_wmv(filename,img_data,x,y,fps,img_func,vtitle)
% makes a movie of a number of datasets after each other
if nargin < 7
    vtitle = '';
end

if nargin < 6
    img_func = @image;
end

if nargin < 5
    fps = 25;   
end

if nargin < 4
    x = [1 size(img_data{1},2)];
end

%aviObj = avifile(filename,'fps',fps,'Quality',100,'compression','Cinepak','Colormap',gray(256));
%aviObj = avifile(filename,'fps',fps,'Quality',100,'compression','None','Colormap',gray(256));

%Movie settings

% p.qualityWmv = 98;
% p.width = 640;
% p.height = 600;
% p.figNo = 100;
% p.outputDir = '';
% p.outputName = '';
% p.sideBySide = 0;
% p.useGcMatSettings = 1;
% p.frameRateMovie = 30;

% number of datasets
data_count = length(img_data)

% Dummy parameters
Tf=1/fps;
timesVectMovie = [];
setidx = [];
for i = 1:data_count
    timesVectMovietmp = 0:Tf:((size(img_data{i},3)-1)*Tf);
    if ~isempty(timesVectMovie)
        timesVectMovie = [timesVectMovie timesVectMovietmp+Tf+timesVectMovie(end)];
    else
        timesVectMovie = [timesVectMovie timesVectMovietmp];        
    end;
    setidx = [setidx ones(1,length(timesVectMovietmp))*i];
end;

v = struct();
v.height = size(img_data{1},1);
v.width = size(img_data{1},2);
v.times = timesVectMovie;    
conf = struct();
confWmv = conf; 
confWmv.videoQuality = 99;
confWmv.videoCompressor = 'WMVideo9 Encoder DMO';

h_fig = figure;
set(h_fig,'DoubleBuffer','on');
ax = axes('Parent',h_fig);
% set(hfig,'visible','off');
set(ax,'Units','Normalized','NextPlot','replace','Visible','on');
set(h_fig,'color',[0 0 0]);

colormap(gray(256))
scrsz = get(0,'ScreenSize');
set(h_fig,'MenuBar','none','Resize','off','Position',[scrsz(3)/2-v.width/2, 2*scrsz(4)/3-2*v.height/3,v.width,v.height]);

lastSet = 0;
for k=1:length(timesVectMovie) 
    if lastSet ~= setidx(k)
        currentSet = img_data{setidx(k)};
        lastSet = setidx(k);
        frameCounter = 1;
    else
        frameCounter = frameCounter +1;
    end;
    
    img_func(x,y,currentSet(:,:,frameCounter));
    set(ax,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1],'YAxisLocation','right')
    axis(ax,'image');
    
    title(ax,vtitle{setidx(k)},'color',[1 1 1],'Interpreter','None')
    
%     T = get(ax,'TightInset');
%     set(ax,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
%     
    F = getframe(h_fig);
    v.frames(k) = F;
    pause(0.05)
end
% close figure;
close(h_fig);
mmwrite(filename,v,confWmv);