function score = CheckEllipse(DATA, C, varargin)
%AllV.CheckEllipse  Calculate how well an ellipse separtes points

if isfield(C,'ellipse')
    xy = DATA.xy{C.cluster};
    [~,score] = FindEllipse(xy,DATA.clst,'cluster',C.cluster+1,'eval',C.ellipse.xy);
    id = find(DATA.clst == C.cluster+1);
end