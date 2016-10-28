function DockFigures(F, varargin)

for j = 1:length(F)
    set(F(j),'windowstyle','docked');
end