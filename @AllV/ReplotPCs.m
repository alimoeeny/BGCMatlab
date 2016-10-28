function DATA = ReplotPCs(DATA,E, varargin)%AllV.ReplotPCs(DATA,[])  Plot scatterplots for current cluster       %AllV.ReplotPCs(DATA,E)  Plot scatterplots using boundary in E%                 ,...,'string',str)  appends str to the title line (middle of plot)DATA = GetDataFromFig(DATA);setids = [];meanpos = [];autospace = 0;forcetofront = 0;addstr = '';args = {};j = 1;while j <= length(varargin)    if strncmpi(varargin{j},'autospace',6)        autospace = 1;    elseif strncmpi(varargin{j},'setid',5)        j = j+1;        setids = varargin{j};    elseif strncmpi(varargin{j},'showpos',7)        j = j+1;        meanpos = varargin{j};        j = j+1;        dimpos = varargin{j};    elseif strncmpi(varargin{j},'string',5)        j = j+1;        addstr = varargin{j};    elseif strncmpi(varargin{j},'tofront',5)        forcetofront = 1;    end    j = j+1;endif DATA.interactive < 0    return;endif autospace    DATA.plottype = AllV.WhichPlotType(DATA.cluster, DATA.currentcluster);endif DATA.plottype == 11    plotspaces(1:8) = 1;else    plotspaces(1:8) = DATA.plottype;endif isempty(E)    cluster = DATA.cluster;else    cluster = E;endif forcetofront    figure(DATA.toplevel);else    AllV.SetFigure(DATA.tag.top, DATA);endclusterplot = [];if isfield(DATA,'dipvals')    dipvals = DATA.dipvals;else    dipvals = zeros(1,8);endif DATA.usebmi == 0dipvals = zeros(1,8);endif DATA.plottype == 2;    vpts = DATA.vpts;elseif strcmp(DATA.plottype,'adcdvdy');    vpts = DATA.dvpts;    DATA.dV = AllV.GetdVdY(DATA);elseif DATA.plottype == 8;    vpts = DATA.dvpts;    AllVoltages=getappdata(DATA.toplevel,'AllVoltages');    DATA.dV= diff(AllVoltages,1,2);elseif DATA.plottype == 9    vpts = [1 1 2 2; 1 1 3 3];    vpts = repmat(vpts,4,1);end%need to get boundary anyway in case there are other clustersif isempty(E) && isfield(DATA,'cluster') && isfield(DATA.cluster,'space')    if DATA.cluster.space(1) == DATA.plottype %cluster is in this space        E = AllV.BoundaryFromCluster(E,DATA.cluster, DATA.currentcluster);    else        E = AllV.BoundaryFromCluster(E,DATA.cluster,DATA.currentcluster);    endendLabels = {};plots = DATA.pcplots(1:8,:);cspace = zeros(size(plots,1),1);if ischar(DATA.plottype)elseif DATA.plottype == 3 || DATA.plottype ==4    if DATA.usestdtemplates        if DATA.plottype == 15            plots = DATA.tmplots(33:40,:);        elseif DATA.plottype == 4            plots = DATA.tmplots(25:32,:);        else            plots = DATA.tmplots(17:24,:);        end    elseif DATA.plottype == 4        plots = DATA.tmplots(9:16,:);    else        plots = DATA.tmplots(1:8,:);    end    if isfield(DATA,'tmpdips')    dipvals = DATA.tmpdips;    endelseif DATA.plottype == 12 %Try new scores    tpt = find(DATA.spts == 0);    p = DATA.probe(1);    AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');    tvals = squeeze(AllVoltages(p,tpt,:));    pcs(:,1) = squeeze(AllVoltages(p,tpt+6,:))-tvals;    pcs(:,2) = squeeze(AllVoltages(p,tpt+19,:))-tvals;    pcs(:,3) = squeeze(AllVoltages(p,tpt-6,:))-tvals;    if p < size(AllVoltages,1)        pcs(:,4) = squeeze(AllVoltages(p+1,tpt,:));    elseif p > 1        pcs(:,4) = squeeze(AllVoltages(p-1,tpt+8,:));    end    if p > 1        pcs(:,5) = squeeze(AllVoltages(p-1,tpt,:));    elseif p< size(AllVoltages,1)        pcs(:,5) = squeeze(AllVoltages(p+1,tpt+8,:));    else        pcs(:,5) = squeeze(AllVoltages(p,tpt+8,:));    end    plots = DATA.pcplots;elseif DATA.plottype == 13 %show best space according to mahalelseif DATA.plottype == 16 %Templates Plot 3        plots = DATA.tmplots(41:48,:);    elseif DATA.plottype == 11 %show spaces needed for all cells    findagain = 0;    pcs = DATA.pcs;    cls = unique(DATA.clst);    cls= cls(cls>1)-1;    plots = DATA.pcplots;    nplots = 0;        if DATA.cluster.space(1) == 1            npc = min([size(DATA.pcs,2) 20]);            Scores = DATA.pcs(:,1:npc);            plots = DATA.pcplots;            plotspaces(1:8) = 1;        elseif DATA.usestdtemplates            plotspaces(1:8) = 3;            Scores = DATA.TemplateScores;            plots = DATA.tmplots(17:24,:);                    else            plotspaces(1:8) = 3;            Scores = DATA.TemplateScores;            plots = DATA.tmplots(1:8,:);        end        cspace = ones(size(plotspaces)) .* DATA.currentcluster;        nx = 0;        xplots = [];    for j = cls(:)'        C = AllV.GetClusterDef(cluster,j);        if isfield(DATA.cluster,'autofiti') && DATA.cluster.autofiti(1) > 1            newspace = 1;        else            newspace = 0;        end        if j > length(DATA.xy) || isempty(DATA.xy{j}) || newspace            [x,y] = AllV.GetClusterXYData(DATA, C.space([2 3]), 'cluster',j);             xy = cat(2,x(:), y(:));        else            xy = DATA.xy{j};        end        plots(j,1) = 2*(j-1)+1;        plots(j,2) = 2*(j-1)+2;        plots(j,1) = C.space(2);        if length(C.space) > 2            plots(j,2) = C.space(3);        else            plots(j,2) = 1;        end        plotspaces(j) = C.space(1);        if DATA.cluster.auto %only one set of template scores            cspace(j) = 1;        else            cspace(j) = j;        end        if plotspaces(j) == 3 && DATA.currentcluster ~= j            plotspaces(j) = 4;        end        isolation = CalcIsolation(xy, DATA.clst, j+1);        if C.space(1) == 3            Labels{j} = sprintf('%s vs %s: %.1f(%.1f)',DATA.TemplateLabels{C.space(2)}, DATA.TemplateLabels{C.space(2)}, isolation([1 3]));        elseif C.space(1) == AllV.USERSPACE            if isfield(C,'Variables') && ~isempty(C.Variables)                Labels{j} = sprintf('%s vs %s',C.Variables{1},C.Variables{2});                plotnames{j} = C.Variables;            else                Labels{j} = '?';            end        elseif C.space(1) == 1            Labels{j} = sprintf('%d vs %d: %.1f(%.1f)',C.space(2),C.space(3), isolation([1 3]));        end        nplots = j;        if findagain            if j == 1                DATA.cluster.ellipse.xy = FindEllipse(pcs(:,1:2),DATA.cluster.clst,'cluster',j);                DATA.cluster.ellipse.space = C.space;            else                DATA.cluster.next{j-1}.ellipse.xy = FindEllipse(pcs(:,end-1:end),DATA.cluster.clst,'cluster',j);            end            cluster = DATA.cluster;        end%if bestiolation space is not the space used to cut, add this too.%but not for manual extra dimensions - need different logic        if isfield(C,'bestisolation') && isfield(C.bestisolation,'space')            if (C.bestisolation.space(1) ~= C.space(1) || sum(~ismember(C.bestisolation.space(2:end),C.space(2:end)))) ...                    && C.bestisolation.space(1) ~= AllV.USERSPACE                nx = nx+1;                xplots(nx,:) = C.bestisolation.space(1:3);                lb = AllV.SpaceLabels(DATA,C,xplots(nx,:));                fprintf('Adding %s vs %s for Cluster %d\n',lb{1},lb{2},j);                if j ~= DATA.currentcluster && xplots(nx,1) ==3                    xplots(nx,1) = 4;                end                xspace(nx) = j;            end        end    end    if nx > 0        plots(j+1:j+nx,:) = xplots(:,2:3);        plotspaces(j+1:j+nx) = xplots(:,1);        cspace(j+1:j+nx) = xspace;        for k = 1:length(nx)            Labels{j+nx} = 'Extra';        end    end    j = j+nx;    dims = unique(plots(1:j,:));%now try and fill up the list with useful plots. Use other pairings%of dimensions already required    while j < 8        for k = 1:length(dims)            for m = 1:k-1                got = sum(ismember(plots(1:j,:),dims([m k])),2);                if max(got) < 2                    j = j+1;                    plots(j,:) = dims([m k]);                end            end        end        j = j+1;    end    if j > 8        plots = plots(1:8,:);    end    if isempty(pcs)        pcs = Scores;    end    npc = size(Scores,2);elseif DATA.plottype == 7    plots = DATA.tmplots(17:24,:);elseif DATA.plottype == 5    plots = DATA.tmplots(9:16,:);elseif DATA.plottype == 1    Labels = AllV.PCLabels(DATA);    elseif DATA.plottype == 14    Labels = AllV.PCLabels(DATA);    elseif DATA.plottype == 0 %% dummy not sure what 0 is    Labels = AllV.PCLabels(DATA);    elseif DATA.plottype == 13    plots = [5 6; 1 2; 1 3; 1 4; 2 3 ; 2 4; 3 4; 1 5;];end   if DATA.usegmcid & length(DATA.gmcid)    if isfield(DATA,'gmcid')        clid = DATA.gmcid;    elseif isfield(E,'bestcl')        clid = E.bestcl;    else    clid = DATA.clid(DATA.uid);    DATA.usegmcid = 0; %don't use it if not defined    endelseif length(setids)    clid = setids;    args = {args{:} 'ptsz' [8 1]};else    args = {args{:} 'ptsz' DATA.ptsz};    if isfield(DATA,'clst') %can be empty if recluste ==4 - don't use uid        clid = DATA.clst;    else        clid = DATA.clid;    endendp = DATA.probelist(DATA.probe(1));type = double(DATA.clplot);if DATA.clplot ==1    if DATA.plot.scaledensity        type = 2;    else        type(2) = DATA.plot.densityscale;    endendC = AllV.GetClusterDef(cluster, DATA.currentcluster);args = {args{:} 'clusterdef' 0};for j = 1:size(plots,1)    mysubplot(2,4,j);    clusterplot(:,j) = AllV.GetClusterPlots(DATA,cluster,plots,j);    if size(clusterplot,1)  < DATA.currentcluster%currentcluster no space defined yet        args{end} = 0;    elseif clusterplot(DATA.currentcluster,j) == j %this is the plot defining the cluster        args{end} = 1;    else        args{end} = 0;    end    t = -1; %not a handle    if DATA.plottype == 11 && plotspaces(j) == 4%when first add a manual template cluser, can have others still using stdtemplates                T = GetTemplateScores(DATA, cspace(j));        AllV.PlotPCs(T.Scores,plots(j,1),plots(j,2),type,clid,DATA.colors,C,'fixrange',args{:});        set(gca,'UserData',[3 plots(j,:)]);        tstr = sprintf('%s vs %s',DATA.TemplateLabels{plots(j,1)},DATA.TemplateLabels{plots(j,2)});        t =  title(tstr);    elseif plotspaces(j) == AllV.USERSPACE        X = AllV.GetValues(DATA, plotnames{j});        AllV.PlotPCs(X,1,2,type,clid,DATA.colors,C,'fixrange',args{:});        t =  title(Labels{j});        clear X;        X.Variables = plotnames{j};        X.ClusterSpace = [AllV.USERSPACE 3 3];        if isfield(C,'space')            X.pcplot = C.space(2:end);        else            X.pcplot = plots(j,:);        end        set(gca,'UserData',X);    elseif ismember(DATA.plottype, [3 4 7 16]) || plotspaces(j) == 3        T = GetTemplateScores(DATA, cspace(j));        AllV.PlotPCs(T.Scores,plots(j,1),plots(j,2),type,clid,DATA.colors,C,'fixrange',args{:});        tstr = sprintf('%s vs %s',DATA.TemplateLabels{plots(j,1)},DATA.TemplateLabels{plots(j,2)});        t =  title(tstr);        if isfield(E,'gmfit') && E.space(1) == 6 && E.space(2) ==4             [a,xi] = ismember(plots(j,1),DATA.tmplspace(1,:));            [b,yi] = ismember(plots(j,2),DATA.tmplspace(1,:));            if a && b && size(E.gmfit.mu,2) > max([xi yi])                hold on;                plot(E.gmfit.mu(1,xi),E.gmfit.mu(1,yi),'g+','markersize',10,'linewidth',2);                plot(E.gmfit.mu(2,xi),E.gmfit.mu(2,yi),'g+','markersize',10,'linewidth',2);            end                    end        set(gca,'UserData',[3 plots(j,:)]);        if isfield(DATA,'oldclusterpts')            hold on;            plot(DATA.TemplateScores(DATA.oldclusterpts{1},plots(j,1)),DATA.TemplateScores(DATA.oldclusterpts{1},plots(j,2)),'+');        end    elseif ismember(DATA.plottype,[6])        AllV.PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),type,clid,DATA.colors);        set(gca,'UserData',vpts(j,:));        if DATA.usegmcid && sum(ismember([vpts(1,3) vpts(3,4)],DATA.vspace) == 2)        end    elseif ismember(DATA.plottype,[2 6 8 9]) | strcmp(DATA.plottype,'adcdvdy')        AllV.PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),type,clid,DATA.colors,'fixrange');        set(gca,'UserData',vpts(j,:));        if DATA.usegmcid && sum(ismember([vpts(j,2) vpts(j,4)],DATA.vspace)) == 2            k = find(DATA.vspace == vpts(j,2));            m = find(DATA.vspace == vpts(j,4));            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)        end    elseif ismember(DATA.plottype,[13]) %plot differeentiated SDC        AllV.PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),type,clid,DATA.colors,'fixrange');        set(gca,'UserData',vpts(j,:));        if DATA.usegmcid && sum(ismember([vpts(j,2) vpts(j,4)],DATA.vspace)) == 2            k = find(DATA.vspace == vpts(j,2));            m = find(DATA.vspace == vpts(j,4));            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)        end    else        if DATA.usegmcid && sum(ismember(DATA.pcplots(j,:),[1:4])) == 2 && isfield(DATA.cluster,'gmfit')            if strcmp(DATA.cluster.autocutmode,'ecker')            showclustermeans = 0;            else            showclustermeans = 1;            end        else            showclustermeans = 0;        end        if ismember(DATA.plottype,[11 12])%tests            if plotspaces(j) == 1                AllV.PlotPCs(pcs,plots(j,1),plots(j,2),type,clid,DATA.colors, C,'fixrange',args{:});            elseif plotspaces(j) == 3                AllV.PlotPCs(DATA.TemplateScoers,plots(j,1),plots(j,2),type,clid,DATA.colors, C,'fixrange',args{:});            end            if length(Labels) >= j && ~isempty(Labels{j})                tstr = Labels{j};            elseif plotspaces(j) == 3                tstr = sprintf('%s vs %s',DATA.TemplateLabels{plots(j,1)},DATA.TemplateLabels{plots(j,2)});            else                tstr = sprintf('%d vs %d',plots(j,1),plots(j,2));            end            t =  title(tstr);        elseif ismember(DATA.plottype,[13]) %james autocut            AllV.PlotFeatures(DATA,plots(j,1),plots(j,2),type,clid,DATA.colors,C, 'fixrange');            showclustermeans = 0;        elseif ismember(DATA.plottype,[10]) %ICA            AllV.PlotPCs(DATA.icas,DATA.pcplots(j,1),DATA.pcplots(j,2),type,clid,DATA.colors, C,'fixrange',args{:});        else            AllV.PlotPCs(DATA.pcs,DATA.pcplots(j,1),DATA.pcplots(j,2),type,clid,DATA.colors, C,'fixrange',args{:});            t = title(Labels{j});        end        set(gca,'UserData',[plotspaces(j) plots(j,:)]);        if showclustermeans            k = find(DATA.pcspace == DATA.pcplots(j,1));            m = find(DATA.pcspace == DATA.pcplots(j,2));            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)        end        if j == 1            text(0,1.2,sprintf('E%dP%d',DATA.exptno,p),'units','normalized','fontsize',DATA.gui.fontsize(1));        end    end    if ~ismember(j,[1 4])        set(gca,'ytick',[]);    end    if ~ismember(j,[5:8])        set(gca,'xtick',[]);    end    if isempty(setids)    axis('tight');              end    xl = get(gca,'Xlim');    yl = get(gca,'Ylim');    if double(t) > 0 && ishandle(t)        set(t,'VerticalAlignment','top','position',[mean(xl) yl(2)]);    end    setappdata(gca,'axiscluster',0);    if DATA.showdipvals    text(mean(xl),yl(2),sprintf('%.2f',dipvals(j)),'VerticalAlignment','top','color','k');    end    %would like to used cluster, not E in future, so it doesn't depend on    %which is current cluster%    clusterplot(:,j) = AllV.GetClusterPlots(DATA,DATA.cluster,plots,j);endif isfield(DATA.Expt,'Header') && isfield(DATA.Expt.Header,'expname')       exname = DATA.Expt.Header.expname;   else       exname = [];   endif DATA.plottype == 1    p = DATA.probelist(DATA.probe);    if isfield(DATA,'alldips')        if size(DATA.alldips,2) > 1            b = max(DATA.alldips');        else            b = DATA.alldips;        end        tstr = sprintf('%s %.2f, dt%.2f, csd%.2f',num2str(p),b(1).*100,b(2).*100,b(3).*100);    elseif DATA.dvdt        tstr = sprintf('E%dP%d: dvdt %s',DATA.exptno,p(1),exname);    elseif DATA.csd        tstr = sprintf('E%dP%d: csd %s',DATA.exptno,p(1),exname);    else        tstr = sprintf('E%dP%d %s',DATA.exptno,p(1),exname);    endelse   if DATA.cluster.auto == 3       xstr = sprintf(' Fit %.1f',AllV.GetValue(DATA.cluster,'fitnumber'));   else       xstr = '';   end   tstr = sprintf('E%dP%d%s %s%s%s',DATA.exptno,p(1),DATA.probelabel,exname,xstr,addstr);endif DATA.triggerset > 0    tstr = sprintf('%s Trigger Set %d',tstr,DATA.triggerset+1);endmysubplot(2,4,1);th = text(0.5,0.0,tstr,'units','normalized','VerticalAlignment','Top','fontsize',DATA.gui.fontsize(1));if isfield(DATA,'cluster') && isfield(DATA.cluster,'auto') && DATA.cluster.auto == 0    set(th,'color','r');endset(th,'Tag','PCTitleString');DATA.maintitle = th;DATA.clustericon = AllV.SetClusterIcon(DATA);    Cs = cluster;for j = 1:size(clusterplot,2)for k = 1:size(clusterplot,1)if clusterplot(k,j)    mysubplot(2,4,clusterplot(k,j)); %need to find right graph    hold on;    set(gca,'xlimmode','manual','ylimmode','manual'); %for 2014b    setappdata(gca,'axiscluster',find(clusterplot(:,j)>0));    if k > 1        Cs = cluster.next{k-1};        Cs.color = DATA.colors{k+1};        if isfield(Cs,'optellipse')            h = DrawEllipse(Cs.optellipse.xy);            set(h,'linestyle','--','color',Cs.color);        end        [DATA.elmousept.h, x] = AllV.DrawEllipse(Cs);    else        Cs = cluster;        Cs.color = DATA.colors{2};        if isfield(Cs,'optellipse')            h = DrawEllipse(Cs.optellipse.xy);            set(h,'linestyle','--','color',Cs.color);        end        if ~isnan(Cs.crit)            [DATA.elmousept.h,x] = AllV.DrawEllipse(Cs);        elseif isfield(Cs,'ellipse')            [DATA.elmousept.h,x] = AllV.DrawEllipse(Cs);        else            DATA.elmousept.h = NaN;            x.markpt = [NaN NaN];        end        DATA.elmousept = CopyFields(DATA.elmousept, x);    end    if ~isempty(DATA.elmousept.h) && myhandle(DATA.elmousept.h)        if isfield(Cs,'isSU') && Cs.isSU(1)            set(DATA.elmousept.h,'linewidth',2);        end    end    if isfield(DATA,'oldclusterpts')        if k ==1            hold on;            plot(DATA.xy{1}(DATA.oldclusterpts{1},1),DATA.xy{1}(DATA.oldclusterpts{1},2),'+');        end    end    if isfield(DATA.elmousept,'h') & ishandle(DATA.elmousept.h)        DATA.elmousept.handles(k) = DATA.elmousept.h;    end    if Cs.shape == 1        tmp = Cs;        tmp.pos(1) = mean(E.pos([1 3])) + diff(E.pos([2 4]))/2;        tmp.pos(3) = mean(E.pos([1 3])) - diff(E.pos([2 4]))/2;        tmp.pos(2) = mean(E.pos([2 4])) - diff(E.pos([1 3]))/2;        tmp.pos(4) = mean(E.pos([2 4])) + diff(E.pos([1 3]))/2;        h = AllV.DrawEllipse(tmp,'r');        set(h,'linestyle',':');    end    hold off; endendendif length(DATA.elmousept.handles) >= DATA.currentcluster;    DATA.elmousept.h = DATA.elmousept.handles(DATA.currentcluster);end[iscell, cellid] =  AllV.isacell(DATA, DATA.exptno, AllV.ProbeNumber(DATA));if iscell    C = DATA.cluster;    mysubplot(2,4,4);    for j = 1:length(cellid)        if cellid(j) > 0        text(1,0.1*j,sprintf('Cell %d',cellid(j)),'units','normalized','color',DATA.colors{j+1},...            'HorizontalAlignment','Right','fontsize',DATA.gui.fontsize(1));        if j > 1 && (length(C.next) < j-1 || isempty(C.next{j-1}))            cprintf('red','Cell Defined for Cluster %d, but the cluster is empty\n',j);        end        end    endendif isfield(DATA,'energy')  && DATA.plot.vare    AllV.SetFigure(DATA.tag.vare, DATA);    subplot(1,1,1);    AllV.PlotVarE(DATA);endif isfield(DATA,'xyplot') && DATA.interactive > 0    if C.space(1) == 8  && isfield(C,'Variables')        AllV.PlotOneXY(DATA, C.Variables,'default');    else        AllV.PlotOneXY(DATA, DATA.xyplot.xy,'default');    endendif nargout == 0    SetData(DATA);endfunction T = GetTemplateScores(DATA, cl)        if DATA.usestdtemplates || cl < 1            T.Scores = DATA.TemplateScores;        else            allT = getappdata(DATA.toplevel,'AllTemplates');            if  cl > length(allT) || isempty(allT(cl).Scores)                cprintf('blue','No Template Scores for Cl%d\n',cl);                T.Scores = DATA.TemplateScores;            else                T = allT(cl);                if DATA.profiling                    fprintf('Templates for Cl%d made at %s\n',cl,T.calctime);                end            end        end