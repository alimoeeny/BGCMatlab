function it = FindChild(F, varargin)
%calls find obj for a window and any paired windows

it = findobj(F,varargin{:});
figs = findobj(allchild(0),'flat','type','figure');
for j = 1:length(figs)
    if isappdata(figs(j),'ParentFigure')
        p = getappdata(figs(j),'ParentFigure');
        if p == F
            a = findobj(figs(j),varargin{:});
            it = [it a];
        end
    end
end
