classdef overlayImgEvalROI < overlayImg
    
    properties
        buttonPressed = false;  % is mouse button pressed?
        paintRoiFct = [];
        roiHdl = [];
    end
    methods
        function obj = overlayImgEvalROI(aBackImg, aOverImg, aOverRange, aTranspRange) 
            obj = obj@overlayImg(aBackImg, aOverImg, aOverRange, aTranspRange);
            set(obj.mFigHdl,'WindowButtonDownFcn',@obj.mouseButtonDownFcn,'Interruptible','off');
            set(obj.mFigHdl,'WindowButtonUpFcn',@obj.mouseButtonUpFcn,'Interruptible','off');
            set(obj.mFigHdl,'WindowButtonMotionFcn',@obj.mouseButtonMotionFcn,'Interruptible','on');
            
            obj.update();
        end
        function mouseButtonMotionFcn(obj,src,evnt) %#ok<INUSD>
            if (obj.buttonPressed && ~isempty(obj.paintRoiFct))
                saveunits = get(obj.mAxHdl,'units');
                set(obj.mAxHdl,'units','pixels');
                axPos = get(obj.mAxHdl,'Position');
                set(obj.mAxHdl,'units',saveunits);
                if strcmp(get(src,'Type'),'figure')
                    % convert coordinates
                    currP = get(src,'CurrentPoint');

                    % check if cursor is inside the axis
                    if ~isInsideRect(currP,axPos)
                        return;
                    end;
                    
                    xLim = get(obj.mAxHdl,'XLim');
                    yLim = get(obj.mAxHdl,'YLim'); 
                    currP = currP - axPos([1 2]); % drop offset
                    currP(2) = axPos(4) - currP(2); % swap y axis
                                
                    currP = currP./axPos([3 4]); % scale to 0..1
                    currP = currP.*[diff(xLim) diff(yLim)];
                    currP = currP + [xLim(1) yLim(1)];
                    if ~isempty(obj.roiHdl) && ishandle(obj.roiHdl)
                        delete(obj.roiHdl);
                    end;
                    hold(obj.mAxHdl,'all');
                    obj.roiHdl = obj.paintRoiFct(currP(1), currP(2),obj.mAxHdl);
                    [oVal, bVal] = obj.giveValues(currP);
                    colIdx = obj.idxImg(cat(1,obj.overImg,repmat(oVal,1,size(obj.overImg,2))), obj.overRange, obj.overMap, obj.overRangeSym);
                    col = obj.overMap(colIdx(end,1),:);                    
                    set(obj.roiHdl,'FaceColor',col,'EdgeColor','w');
                    hold(obj.mAxHdl,'off');
                end;
            end;
        end;
        function mouseButtonDownFcn(obj,src,evnt)
            obj.buttonPressed = true;
            obj.mouseButtonMotionFcn(src,evnt);
%             disp('button down');
        end;
        function mouseButtonUpFcn(obj,src,evnt)
            % forward request to update position
            obj.mouseButtonMotionFcn(src,evnt);
            obj.buttonPressed = false;
%             disp('button up');
            if ~isempty(obj.roiHdl) && ishandle(obj.roiHdl)
                delete(obj.roiHdl);
            end;
        end;
        function setPaintRoiFct(obj, aFcnHdl)
            obj.paintRoiFct = aFcnHdl;
        end;
    end
end