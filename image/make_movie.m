function make_movie(filename,img_data,x,y,fps,img_func,type)

if nargin < 8
    vtitle = '';
end

if nargin < 7
    img_func = @image;
end

if nargin < 6
    img_func = @image;
end

if nargin < 5
    fps = 25;   
end

if nargin < 4
    y = [1 size(img_data,1)];
end

if nargin < 3
    x = [1 size(img_data,2)];
end

% [dummy1,dummy2,ext] = fileparts(filename);
% if ~strcmp(ext,['.' type])
%     filename = [filename,'.',type];
% end

if strcmp(type,'avi')
    aviObj = avifile(filename,'fps',fps,'Quality',100,'Colormap',gray(256),'compression','None');
end

h_fig = figure;
set(h_fig,'Color',[1 1 1]);
units=get(h_fig,'units');
%set(h_fig,'units','normalized','outerposition',[0.1 0.305  1.4*[0.3298  0.4695]]);
get(h_fig,'outerposition');
set(h_fig,'units',units);

set(h_fig,'DoubleBuffer','on');
ax = axes('Parent',h_fig,'visible','off');
set(ax,'Units','Normalized','NextPlot','replace','Visible','on');

colormap(gray(256));

if strcmp(type,'mpeg')
    mpgrate=25;             %MPG has 25 frames per second as default
    rate=mpgrate/fps;       %Calculates updated framerate for obtaining 25 fps
    times=1:1/rate:size(img_data,3);            
else
    times = 1:size(img_data,3);
end

for k=1:length(times)   
    set(get(ax,'parent'),'renderer','zbuffer');
    img_func(x,y,img_data(:,:,round(times(k))));   
    set(ax,'YAxisLocation','right')
    axis(ax,'image');
    axis image;shading interp;
    set(h_fig,'Color',[0 0 0])
    set(gca,'YColor',[1 1 1])
    set(gca,'XColor',[1 1 1])
    h1=get(gca,'Title');
    set(h1,'Color',[1 1 1])
    F = getframe(h_fig);
    
    if strcmp(type,'avi')
        aviObj = addframe(aviObj,F);     
    else
        mov(k) = F;
    end
end

if strcmp(type,'avi')
    aviObj = close(aviObj);
else
    map=[0:255;0:255;0:255]'./255;
    mpgoptions=[1, 0, 1, 0, 10, 5, 5, 5]; % Slow, but good quality
    %  mpgoptions=[1, 1, 1, 0, 10, 8, 10, 25];
    mpgwrite(mov,map,filename,mpgoptions);
end