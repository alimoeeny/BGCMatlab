function M = CellToMat(C, f, varargin)
%CellToMat M = CellToMat(C, f,...
%takes the field f from each element of C, and puts into
% a matrix. Only works for cell arrays with 1 or 2 dimensions
M = [];
if size(C,1) == 1 && size(C,2) > 1
    for k = 1:length(C)
        if iscell(C{k})
            for j= 1:length(C{k})
                if isfield(C{k}{j},f)
                    M(k,j,1:length(C{k}{j}.(f))) = C{k}{j}.(f);
                else
                    M(k,j,:) = NaN;
                end
            end
        elseif isfield(C{k},f)
            M(k,1:length(C{k}.(f))) = C{k}.(f);
        end
    end
else
for j = 1:size(C,1)
    for k = 1:size(C,2)
        if isfield(C{j,k},f)
            M(j,k,:) = C{j,k}.(f);
        end
    end
end
end