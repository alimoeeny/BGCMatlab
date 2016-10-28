function PlotAllCellXY(DATA, varargin)listtype = 'default';j = 1;while j <= length(varargin)    if strcmp(varargin{j},'autolist')        listtype = 'autolist';    end    j = j+1;endmuCellList = [];    cellid = PC.GetCurrentCell(DATA);    if strcmp(listtype,'autolist');        Clusters = getappdata(DATA.toplevel,'AutoClusters');        if isappdata(DATA.toplevel,'AutoCellList')            A = getappdata(DATA.toplevel,'AutoCellList')            CellList = A.CellList;        else            CellList = DATA.autolist.CellList;        end        cells = PC.SelectedCells(DATA.selectautoprobe, CellList);        cellid = cells;        if isfield(DATA.autolist,'muCellList')            muCellList = DATA.autolist.muCellList;        end    else        Clusters = getappdata(DATA.toplevel,'Clusters');        CellList = DATA.CellList;        if isfield(DATA,'muCellList')            muCellList = DATA.muCellList;        end       cellid = PC.SelectedCells(DATA.selectprobe, CellList);    end        ts = now;    if DATA.profiling        profile on;    end    id = find(CellList == cellid);    [eid, cid, clid] = ind2sub(size(CellList),id);    if ~isempty(muCellList)        muid = find(muCellList == cellid);        [mueid, mucid, muclid] = ind2sub(size(muCellList),muid);        eid = cat(1,eid, mueid);        cid = cat(1, cid, mucid);        clid = cat(1, clid, muclid);        allid = cat(1,id, muid);    else        muid = [];        allid = id;     end        [eid, id] = sort(eid);    allid = allid(id);    cid = cid(id);    clid = clid(id);    cid = 1+mod(cid-1,DATA.nprobes);    [nr,nc] = Nsubplots(length(eid));    PC.SetFigure(DATA,DATA.tag.allxy);    subplot(1,1,1);    colors = mycolors('spkcolors');    for j = 1:length(eid)        mysubplot(nr,nc,j);        hold off;         C = PC.GetClusterInfo(Clusters, [eid(j) cid(j) 1]);        plots = PC.PlotClusterXY(DATA,C,'shorttitle','cellid',clid(j),'allxy');        if clid(j) ~= 1 && length(plots) > 1            np = clid(j)+1;            if length(plots) >= np                set(plots(np),'color','r');                set(plots(2),'color',colors{np});            else                set(plots(2),'color','k');            end        end        if DATA.plot.trighist            PC.AddTrigHist(DATA,C, DATA.currentcluster);        end        set(gca,'Xtick',[],'Ytick',[]);        h = get(gca,'title');        xl = get(gca,'xlim');        yl = get(gca,'ylim');        a = get(h,'position');        a(2) = yl(2);        a(1) = mean(xl);        set(h,'position',a,'VerticalAlignment','top');        color = 'b';        if max(C.mahal([1 3])) < 3            color = 'r';        elseif max(C.mahal([1 3])) < 2            color = 'k';        end                   if length(plots) > 1            ph = plots(2:end);        else            ph = [];        end        if ismember(allid(j),muid)            set(h,'fontweight','normal','fontsize',10,'color','k');        else            set(h,'fontweight','normal','fontsize',12,'color',color);        end        set(gca,'ButtonDownFcn',{@PC.HitXYPlot, eid(j), cid(j)});        c = get(gca,'Children');        for k = 1:length(c)            if ~ismember(c(k),ph);            set(c(k),'ButtonDownFcn',{@PC.HitXYPlot, eid(j), cid(j)});            end        end        PC.AddCellLabels(DATA, eid(j), cid(j));    end    if DATA.profiling        fprintf('Took %.2f\n',mytoc(ts));        profile viewer;    end    