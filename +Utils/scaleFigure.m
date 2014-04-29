function scaleFigure(varargin)
    % scales a figure
    % scaleFigure(scaleXY) % scales gcf in both x and y with given factor
    % scaleFigure(scaleX, scaleY) % scales gcf in both x and y with given
    %                   factor
    % scaleFigure(scaleXY,[],figH) % scales figH in both x and y
    % scaleFigure(scaleX, scaleY, figH)
    
    % Jochen, 2010
    
    
    switch nargin
        case 1   % scaleFigure(scaleXY) % scales gcf in both x and y with given factor
            figH = gcf;
            scaleX = varargin{1};
            scaleY = varargin{1};
        case 2   % scaleFigure(scaleX, scaleY) 
            figH = gcf;
            scaleX = varargin{1};
            scaleY = varargin{2};
        case 3    
            if isempty(varargin{2}) % scaleFigure(scaleXY, [], figH)            
                figH = varargin{3};
                scaleX = varargin{1};
                scaleY = varargin{1};
            else                  % scaleFigure(scaleX, scaleY, figH)
                figH = varargin{3};
                scaleX = varargin{1};
                scaleY = varargin{2};
            end
    end
    windecory = 70;
    
    % get current size
    pos = get(figH,'Position');
    
    width = pos(3)*scaleX;
    height = pos(4)*scaleY;
    
    screensize = get(0,'ScreenSize');
    height = min(height, screensize(4)-windecory-5);
    width = min(width, screensize(3)-10);
    
    newpos = floor([pos(1)+(pos(3)-width)/2 pos(2)+(pos(4)-height)/2 width,height]);
    newpos(1) = max(newpos(1),1);
    newpos(2) = min(newpos(2),screensize(4)-height-windecory);
    newpos(2) = max(newpos(2),5);
    
    set(figH,'Position',newpos);
end
            
        