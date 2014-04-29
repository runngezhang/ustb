function make_movie2_wmv(filename,img_data1,img_data2,x,y,fps,img_func,vtitle)

if nargin < 8
    vtitle = '';
end

if nargin < 7
    img_func = @image;
end

if nargin < 6
    fps = 25;   
end

if nargin < 5
    y = [1 size(img_data1,1)];
end

if nargin < 4
    x = [1 size(img_data1,2)];
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

% Dummy parameters
Tf=1/fps;
timesVectMovie = 0:Tf:((size(img_data1,3)-1)*Tf);

v = struct();
% v.height = y;
% v.width = x;
v.height = size(img_data1,1);
v.width = size(img_data1,2)*2;
v.times = timesVectMovie;    
conf = struct();
confWmv = conf; 
confWmv.videoQuality = 99;
confWmv.videoCompressor = 'WMVideo9 Encoder DMO';


h_fig = figure;
set(h_fig,'DoubleBuffer','on');
ax = axes('Parent',h_fig,'visible','off');
set(ax,'Units','Normalized','NextPlot','replace','Visible','on');
set(h_fig,'color',[0 0 0]);

colormap(gray(256))
scrsz = get(0,'ScreenSize');
set(h_fig,'MenuBar','none','Resize','off','Position',[scrsz(3)/2-v.width/2, scrsz(4)/2-v.height/2,v.width,v.height]);
% set(h_fig,'Visible','off');
% set(ax,'Visible','off');
for k=1:size(img_data1,3) 
    img_func(x,y,img_data1(:,:,k),'Parent',ax);
    hold on;  
    img_func(x+diff(x)*1.05,y,img_data2(:,:,k),'Parent',ax);    

    set(ax,'xtick',[0:10:diff(x),(0:10:diff(x))+diff(x)*1.05], ...
        'xticklabel',[0:10:diff(x),0:10:diff(x)]);
    set(ax,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1],'YAxisLocation','right')
    
    %line([diff(x),diff(x)]*1.025,y,'color',[1,1,1])
    
    axis(ax,'image');
    
    title(ax,vtitle,'color',[1 1 1],'Interpreter','None')
    
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