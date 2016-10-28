function [cells,e,p,clid] = SelectedCells(selectprobe, CellList, varargin);
%cells = SelectedCells(selectid, CellList, varargin);
%cells = SelectedCells(DATA);
%generate list of current selected cells
listtype = 'default';

if isstruct(selectprobe)
    DATA = selectprobe;
    listtype = PC.GetValue(DATA,'listtype');
    if strcmp(listtype,'autolist')
        if ~isfield(DATA,'selectautoprobe')
            DATA.selectautoprobe = DATA.selectprobe;
        end
        selectprobe = DATA.selectautoprobe;
        if isfield(DATA.autolist,'CellList')
            CellList = DATA.autolist.CellList;
        else
            CellList = [];
        end
    else
        selectprobe = DATA.selectprobe;
        CellList = DATA.CellList;
    end
end
       [e,p] = find(selectprobe > 0);
        id = sub2ind(size(selectprobe),e,p);
        clid = selectprobe(id);
        cid = sub2ind(size(CellList),e,p,selectprobe(id));
        cells = unique(CellList(cid));
        cells = cells(cells>0);
        if isempty(cells) %may have selected duplicate
            cells = unique(CellList(cid));
            cells = cells(abs(cells) > 0);
        end
end