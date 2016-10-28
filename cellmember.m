function x = cellmember(a,b, varargin)
%cellmember(a,b) ismember for cell arrays of doubles.
%if a is a scalar, b a cell array, returns which cellarrays contain a
%if a is a vector, returns a matrix length(a) x length(b) with
% membership for each element of a
if isnumeric(a) && iscell(b)
    for j = 1:length(b)
        if isnumeric(b{j})
            x(:,j) = ismember(a,b{j});
        end
    end
end