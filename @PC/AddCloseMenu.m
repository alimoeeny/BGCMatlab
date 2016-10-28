function sm = AddCloseMenu(F, DATA, varargin)
%sm - PC.AddCloseMenu(F, DATA, varargin) Add 'close plotclusters' to a menu
sm = [];

fm = findobj(allchild(F),'flat','type','uimenu','tag','figMenuFile');
if ~isempty(fm)
    sm = uimenu(fm,'Label','Close PlotClusters','callback', {@PC.PlotClusters, 'close'});
end

%for generic case....
function closefigure(a,b,F)
close(F);