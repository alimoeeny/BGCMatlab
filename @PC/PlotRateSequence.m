function Expts = PlotRateSequence(DATA, pos, varargin)
%PC.PlotRateSequence(DATA, pos)  plot rate sequence for a series of expts

if isfigure(DATA)
    F = DATA;
    DATA = GetDataFromFig(F);
end

Expts = getappdata(DATA.toplevel, 'Expts');
CellList = PC.GetValue(DATA,'CellList');
C = PC.GetClusterInfo(DATA, pos,'clst');
exid = pos(:,1);
for k = 1:length(C)
    [xcl, b] = PC.FindExcludedTrials(DATA,pos(k,1),pos(k,2),pos(k,3),C{k});
    aExpts{k} = PC.CountSpikes(Expts{pos(k,1)},C{k},pos(k,3),b.xid);    
end
colors = mycolors;
scaling = 'rates';
cells = celllist.find(CellList, pos, 'lookup');
cells = cells(cells>0);
a = PlotRateSequence(aExpts,'color',colors{1},scaling,'bytime','callback',{@PC.HitTrial, cells(1)});
yl = get(gca,'ylim');
gaps = 0;
for j = 1:length(a.AllBlocks)
    text(a.AllBlocks(j),yl(2),sprintf('%d:%s',GetExptNumber(aExpts{j}),Expt2Name(aExpts{j})),'rotation',90,'horizontalalignment','right','verticalalignment','top');
    if j > 1
        gaps(j) = ((Expts{exid(j)}.Trials(1).TrialStart -Expts{exid(j-1)}.Trials(end).TrialStart)./10000)...
            +Expts{exid(j)}.Header.timeoffset-Expts{exid(j-1)}.Header.timeoffset -2;
    end
    Expts{exid(j)}.Header.timeadjust = sum(gaps);
    X.expstarts(j) = Expts{exid(j)}.Header.timeoffset-Expts{exid(j)}.Header.timeadjust ...
        +Expts{exid(j)}.Trials(1).TrialStart/10000;
end
    

X.yrange = get(gca,'ylim');
X.xrange = get(gca,'xlim');
X.exptlist = pos(:,1);
X.cellids = cells;
X.timeadjust = cumsum(gaps);

set(gcf,'UserData',X);