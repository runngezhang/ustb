function [lineHdl] = addVertLine(pos,ax,varargin)
    % addVertLine(pos,ax,varargin)
    % adds a vertical line over the height of the plot
    % pos - position at the x axis where the line will be added
    % ax - axes handle if not given gca is used
    % varargin - further line properties which will be passed to the line
    %            function (color, line-style ...)
    if ~exist('ax','var')
        ax = gca;
    elseif isempty(ax)
        ax = gca;
    end;
    ytmp = get(ax,'ylim');
    lineHdl = line(pos*[1 1],ytmp,'Parent',ax,varargin{:});
end