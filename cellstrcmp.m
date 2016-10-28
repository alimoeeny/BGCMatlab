function [matches, ids] = cellstrcmp(a,b)
%[nmatches, ids] = cellstrcmp(a,b) convenience strcmp with cell array

if iscell(b) && ~iscellstr(b)
    b = cell2cellstr(b);
end
if iscell(a)
    for j = 1:length(a)
        [matches(j), ids(j,:)] = cellstrcmp(a{j},b);
    end
else
    ids = strcmp(a,b);
    matches = sum(ids);
end