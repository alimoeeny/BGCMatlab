function DATA = SetCellEntry(DATA, C,  e, p, c, cellid, varargin)        savelist = 1;    muonly = 0;    setduplicate = 0;    j = 1;     while j <= length(varargin)        if strncmpi(varargin{j},'duplicate',6)            setduplicate = 1;            if length(varargin) > j && isnumeric(varargin{1+1})                setduplicate = varargin{j+1};                j = j+1;            end        elseif strncmpi(varargin{j},'muonly',6)            muonly = 1;        elseif strncmpi(varargin{j},'nosave',6)            savelist = 0;        end        j = j+1;    end    if isempty(C)        Clusters = getappdata(DATA.toplevel,'Clusters');        C = Clusters{e}{p};    end    quality = 4;    if c < 0        id = find(DATA.CellList(e,p,:) == cellid);        if length(id) == 1            quality = 4 + c;            c = id;        else            return;        end    end    if c <= size(DATA.CellList,3)        oldcell = DATA.CellList(e,p,c);    else        oldcell = 0;    end    if cellid == 0 && setduplicate        cellid = DATA.currentcell;    end    if cellid > 0        [a,b] = find(squeeze(DATA.CellList(e,:,:)) == cellid);        if setduplicate            DATA.CellList(e,a,b) = -setduplicate; %remove this cell from other clusters this expt        else            DATA.CellList(e,a,b) = 0; %remove this cell from other clusters this expt        end    end    if muonly        DATA.muCellList(e,p,c) = cellid;    else        if setduplicate == 0            DATA.CellList(e,p,c) = cellid;        end        DATA.muCellList(e,p,c) = 0; %remove from mucelllist    end    DATA.CellChanges = cat(1,DATA.CellChanges,[e p cellid c now oldcell]);    if c > 1        if c > length(C.next)+1            return;        end        C = C.next{c-1};    end    if isfield(C,'excludetrialids') && length(C.excludetrialids)        % find trial #s that match the ids excluded        DATA.CellDetails.excludetrials{e,p,c} = union(DATA.CellDetails.excludetrials{e,p,c}, C.excludetrialids);    end    DATA.CellDetails.Quality(e,p,c) = quality; %default    if savelist        DATA = PC.SaveCellList(DATA);    end