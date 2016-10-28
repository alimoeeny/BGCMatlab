function exid = GetExptid(DATA, e, type, varargin)
%eid = GetExptid(DATA, e, type)
%return index for expt in different listst
%eid{1} is the row number for the CellImage plot
%eid{2} is the element number is the Expts  cell array (DATA.exptlist)
%eid{3} is the element number in the loaded Clusters Array (DATA.exptids)
%
%DATA.exptid is a list of expt numbers in teh Expts Struct
if nargin < 3
    type = 'row';
end
if strcmp(type,'row')
    if isfield(DATA.CellDetails,'exptids')
        eid = DATA.CellDetails.exptids(e);
    else
        eid = DATA.exptid(e);
    end
elseif strcmp(type,'autolist')
    exid{1} = DATA.exptid(e);
    exid{3} = e;
    return;
elseif strcmp(type,'expt')
    eid = e;
elseif strcmp(type,'loaded')
    eid = DATA.exptids(e);    
else
    eid = DATA.CellDetails.exptids(e);
end


%eid is now absolute expt number
x = find(ismember(DATA.CellDetails.exptids,eid));
xid = 1;
if isempty(x)
    exid{1} = NaN;
elseif length(x) == 1
    exid{1} = x;
else
    if DATA.CellDetails.exptids(e-1) == eid
        exid{1} = x(2);
        xid = 2;
    else
        exid{1} = x(1);
    end
end

x = find(ismember(DATA.exptlist,eid));
if isempty(x)
    exid{2} = NaN;
else
    exid{2} = x(xid);
end

x = find(ismember(DATA.exptid,eid));
if isempty(x)
    exid{3} = NaN;
else
    exid{3} = x(xid);
end
