    function PlotCellRates(DATA,type)    %PC.PlotCellrates(DATA,type) plots rate sequence for currentcell    %        Expts = getappdata(DATA.toplevel,'Expts');        Clusters = getappdata(DATA.toplevel,'Clusters');        [e, p, c] = PC.FindCell(DATA, DATA.currentcell);        id = e;        PC.SetFigure(DATA,DATA.tag.xyseq);        for j = 1:length(id)            xcl = PC.FindExcludedTrials(DATA,e(j),p(j),c(j),Clusters{e(j)}{p(j)});            xid = [Expts{id(j)}.Trials(xcl).id];            Expts{id(j)} = PC.CountSpikes(Expts{id(j)},Clusters{e(j)}{p(j)},c(j),xid);        end        plottype = strmatch(type, {'both' 'rates' 'xy'});        if ismember(plottype,[1 2])            if plottype == 1                subplot(2,1,1);                hold off;            else                subplot(1,1,1);                hold off;            end            PlotRateSequence(Expts(id));        end        if ismember(plottype,[1 3])            if plottype == 1                subplot(2,1,2);            else                subplot(1,1,1);            end            hold off;            for j = 1:length(id)                PC.PlotXYSequence(DATA, Clusters{e(j)}{p(j)}, 'expt', Expts{id(j)});                hold on;            end        end        