function [cellid, it] = GetCurrentCell(DATA,varargin)
%cellid = PCGetCurrentCell(DATA,varargin)
% read cell number from gui menu if there, or use DATA.curentcell
cellid = [];
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'set')
        j = j+1;
        cellid = varargin{j};
    end
    j = j+1;
end
if isfield(DATA.fig,'celllist') && isfigure(DATA.fig.celllist)
    F = DATA.fig.celllist;
else
    F = GetFigure(DATA.tag.celllist);
end
it = findobj(allchild(F),'flat','Tag','CellNumberId');
if cellid > 0 %setting
    set(it,'value',cellid);
else
    cellid = get(it,'value');
    if isempty(cellid)
        cellid = DATA.currentcell;
    end
end