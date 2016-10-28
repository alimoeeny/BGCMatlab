function [s, keep] = cell2cellstr(C, varargin)
%s = cell2cellstr(C) make cellstr from cell that may have none-cell elements
%[s, sid] = cell2cellstr(C) returns indices in C that are strings
%s = cell2cellstr(C,field)  make cell string array using field from each
%cell 
s = {};
if ~iscell(C) ||  isempty(C)
    return;
end
j = 1;
while j <= length(varargin)
    if ischar(varargin{j})
        [s, keep] = GetStrFromCell(C, varargin{j});
        return,
    end
    j = j+1;
end

keep = [];
for j = 1:length(C)
    if ischar(C{j})
        keep(j) = 1;
    end
end
s = C(keep>0);


function [s, keep] = GetStrFromCell(C, f)

s = {};
for j = 1:length(C)
    if isfield(C{j},f)
        s{j} = C{j}.(f);
        keep(j) = 1;
    else
        keep(j) = 0;
    end
end
s = s(keep(keep>0));