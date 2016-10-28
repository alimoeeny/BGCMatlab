function [ok, matches] = iscellfield(X, f, vararagin)
%iscellfield(X, f, vararagin) does cell array have structures with field X
ok =0;
matches = 0;
if iscell(X)
    for j = 1:length(X)
        if isfield(X{j},f)
            matches(j) = 1;
        else
            matches(j) = 0;
        end
    end
ok = sum(matches) > 0;
end
