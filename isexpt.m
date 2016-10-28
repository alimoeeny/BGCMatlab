function yes = isexpt(E, varargin)
%isexpt(E) true if E is an Expt Struct, or cell array
if iscell(E)
    for j = 1:length(E)
        yes(j) = isexpt(E{j});
    end
elseif isfield(E,'Header') && isfield(E,'Trials')
    yes = true;
else
    yes = false;
end