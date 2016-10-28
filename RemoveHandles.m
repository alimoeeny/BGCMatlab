function X = RemoveHandles(X, varargin)
%RemoveHandles(X) removes figure handles from structure X
%so that can save without them


if isstruct(X) && length(X) == 1
    f = fields(X);
    for j = 1:length(f)
        if isfigure(X.(f{j}))
            X = rmfield(X,f{j});
        elseif isstruct(X.(f{j}))
            X.(f{j}) = RemoveHandles(X.(f{j}));
        end
    end
elseif iscell(X)
    for j = 1:length(X)
        X{j} = RemoveHandles(X{j});
    end
end
