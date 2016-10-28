function varargout = SetValue(DATA, type, value, varargin)
%PC.SetValue(DATA, type, ...)  Sets properties safely.
%PC.SetValue(DATA, 'Expt', ex) Sets DATA.Expt to expt matcing row ex
%PC.SetValue(DATA, 'nclusters') recalculates DATA.nclusters
DATA = GetDataFromFig(DATA);

if strcmp(type,'Expt')
    Expts = getappdata(DATA.toplevel,'Expts');
    eid = floor(DATA.exptlist(value));
    DATA.currentpoint(1) = value;
    if length(Expts) >= eid
        DATA.Expt  = Expts{eid};
    end
    varargout{1} = DATA;
elseif strcmp(type,'nclusters')
    C = PC.GetValue(DATA, 'Clusters');
    for j = 1:length(C)
        for k = 1:length(C{j})
            DATA.nclusters(j,k) = max(C{j}(k).clst);
        end
    end
elseif strcmpi(type,'celllist')
    listtype = PC.GetValue(DATA,'listtype');
     if strcmp(listtype,'autolist')
         DATA.autolist.CellList = value;
     else
         DATA.CellList = value;
     end
elseif strcmpi(type,'oldcelllist') %made before any fixes
    C = getappdata(DATA.toplevel,'AutoCellListBackup');
    DATA.autolist.CellList = C;
elseif strcmpi(type,'lastcelllist') %made before each fix
    C = getappdata(DATA.toplevel,'LastCellList');
    if ~isempty(C)
    DATA.autolist.CellList = C;
    end
end


if nargout == 0
    SetData(DATA);
else
    varargout{1} = DATA;
end
