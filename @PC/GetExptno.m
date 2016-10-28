function exid = GetExptno(DATA, e, type, varargin)
%eid = GetExptno(DATA, eid, type)
%return exptno for an exot indexor expt in different listst
%eid{1} is the row number for the CellImage plot
%eid{2} is the element number is the Expts  cell array (DATA.exptlist)
%eid{3} is the element number in the loaded Clusters Array (DATA.exptids)

if nargin < 3
    type = 'row';
end

if strcmp(type,'row')
    if isfield(DATA.CellDetails,'exptids')
        exid = DATA.CellDetails.exptids(e);
    else
        exid = DATA.exptid(e);
    end
elseif strcmp(type,'expt')
    exid = e;
elseif strcmp(type,'loaded')
    exid = DATA.exptids(e);    
else
    exid = DATA.CellDetails.exptids(e);
end
return;

%eid is now absolute expt number
x = find(ismember(DATA.CellDetails.exptids,eid));
if isempty(x)
    exid{1} = NaN;
else
    exid{1} = x;
end

x = find(ismember(DATA.exptlist,eid));
if isempty(x)
    exid{2} = NaN;
else
    exid{2} = x;
end

x = find(ismember(DATA.exptid,eid));
if isempty(x)
    exid{3} = NaN;
else
    exid{3} = x;
end
