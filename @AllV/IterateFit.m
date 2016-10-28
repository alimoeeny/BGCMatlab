function DATA = IterateFit(DATA, niter)j = 1;xc = [0 0];if DATA.currentcluster == 1    C = DATA.cluster;else    C = DATA.cluster.next{DATA.currentcluster-1};end  AllV.ReplotPCs(DATA,[]);for j = 1:niter    DATA.MeanSpike = AllV.PlotMeanSpike(DATA,'recalc');    AllV.TemplatePlot(DATA,'nodip','noplot');    DATA = get(DATA.toplevel,'UserData');    T{j} = DATA.MeanSpike.ms;    if C.space(1) == 6        G = DATA.cluster.gmfit;%        xy = AllV.ProjectND(DATA, size(G.mu,2), G);        AllV.OptimizeBoundary(DATA);    elseif C.space(1) == 3        DATA = AllV.OptimizeBoundary(DATA);        AllV.ReplotPCs(DATA,[]);    endendset(DATA.toplevel,'UserData',DATA);    