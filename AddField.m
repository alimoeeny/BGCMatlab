function S = AddField(S, name, default, varargin)
% S = AddField(S, name, default) add field name to struct S, if it does not exist
% if default is not given, 0 is used

if iscellstr(name)
    if nargin > 2
        varargin = {default varargin{:}};
    end
    for j = 1:length(name)
        S = AddField(S,name{j}, varargin{:});
    end
    return;
end
if isfield(S,name)
    return;
else 
    if nargin > 2
        S.(name) = default;
    else
        S.(name) = 0;
    end
end