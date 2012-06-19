function [count, vals] = Counts(x, varargin)
%[counts, vals] = Counts(x) 
%return counts for each unique value of x
%[counts, vals] = Counts(x, vals)
%return counts for each value in vals. N.B. not a histogram, 
% only counts exact matches


if isempty(x)
    count = [];
    vals  = NaN;
    return;
end
j = 1;
if length(varargin) & isnumeric(varargin{1})
    vals = varargin{1};
    j = j+1;
else
vals = unique(x);
end

for j = 1:length(vals)
    count(j) = sum(x(:) == vals(j));
end