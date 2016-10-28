function SetParent(F, P, varargin)

if isfigure(P)
    setappdata(F,'ParentFigure',P);
elseif isfield(P,'toplevel')
    setappdata(F,'ParentFigure',P.toplevel);
end
    