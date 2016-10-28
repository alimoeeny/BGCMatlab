function DATA = CheckCellList(DATA,eid,p, varargin)
%PC.CheckCellList(DATA,eid,p, varargin)

%make sure CellList,muCellList dimensinos are up to date
if isfield(DATA,'muCellList')
    sz = size(DATA.CellList);
    if size(DATA.muCellList,3) < sz(3) || size(DATA.muCellList,1) < sz(1) || size(DATA.muCellList,2) < sz(2) 
        DATA.muCellList(sz(1),sz(2),sz(3)) = 0;
    end
    np = max([p DATA.nprobes]);
    if size(DATA.muCellList,2) < np
        DATA.muCellList(end,np,end) = 0;
    end
    if size(DATA.muCellList,1) < eid
        DATA.muCellList(eid,np,end) = 0;
    end
end
if ~isfield(DATA,'nclusters')
    DATA.nclusters = 1;
end
if ~isfield(DATA,'dropi')
    DATA.dropi = 0;
end
if size(DATA.nclusters,1) < size(DATA.CellList,1)
    for j = size(DATA.nclusters,1):size(DATA.CellList,1);
        for p = 1:size(DATA.nclusters,2)
            DATA.nclusters(j,p) = sum(squeeze(DATA.CellList(j,p,:)) > 0);
        end
    end
end

if size(DATA.dropi,1) < size(DATA.CellList,1)
    DATA.dropi(size(DATA.CellList,1),1) = 0;
end
if isfield(DATA.CellDetails,'exptids')
    id = setxor(DATA.exptid,DATA.CellDetails.exptids);
    if ~isempty(id)
        cprintf('red','Mismatch between Expt Lists %s\n',sprintf(' %d',id));
    end
end