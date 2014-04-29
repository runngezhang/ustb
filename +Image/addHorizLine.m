function [lineHdl] = addHorizLine(pos,ax,varargin)
    % addHorizLine(pos,ax,varargin)
    % adds a horizontal line over the width of the plot
    % pos - position at the y axis where the line will be added
    % ax - axes handle if not given gca is used
    % varargin - further line properties which will be passed to the line
    %            function (color, line-style ...)
    if ~exist('ax','var')
        ax = gca;
    elseif isempty(ax)
        ax = gca;
    end;
    xtmp = get(ax,'xlim');
    lineHdl = line(xtmp,pos*[1 1],'Parent',ax,varargin{:});
end