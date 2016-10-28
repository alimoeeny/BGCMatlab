function Cut = PlotHistogram(DATA, E, varargin)Cut.gmdprime = NaN;        plotdip = 0;    checkdprimes = 0;    plotgm = 0;   replotpts = 0;   quickmode.quick = 0;    j = 1;    while j <= length(varargin)        if isfield(varargin{j}, 'quickest')            quickmode = varargin{j};            quickmode.quick = 1;        elseif strncmpi(varargin{j},'plotdip',5)            plotdip = 1;        elseif strncmpi(varargin{j},'plotgmdetails',10)            plotgm = 2;        elseif strncmpi(varargin{j},'plotgm',5)            plotgm = 1;        elseif strncmpi(varargin{j},'quick',5)            quickmode.quick = 1;            plotgm = 0;        elseif strncmpi(varargin{j},'fit1cut',6)            quick = 2;            plotgm = 0;        end        j = j+1;    end% E.shape == 2 means that DATA.xy{1} has already been rotated, so don't rotate here.% if classify spike hasn't been called yet, this doesn't work.  Need to% check for this....if isempty(E)    E = AllV.BoundaryFromCluster([],DATA.cluster, DATA.currentcluster);endif DATA.currentcluster == 1    C = DATA.cluster;else    C = DATA.cluster.next{DATA.currentcluster-1};endif isfield(E,'pcplot')    pcplot = E.pcplot;elseif isfield(E,'space')    pcplot = E.space(2:end);endif ~isfield(E,'shape') %can happen if onlye GM defined    return;endif E.shape(1) == 2 && strcmp(E.autocutmode,'ecker')    return;endif E.shape(1) == 3    return;endif E.shape(1) == 2 || E.space(1) == 6%    xy = DATA.xy{1};    angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));    exy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);    crit = mean(exy(:,1));    xy = xyrotate(DATA.ndxy(:,1),DATA.ndxy(:,2),angle);else    pos = E.pos;[allx, ally] = AllV.GetClusterXYData(DATA, pcplot);    if E.shape == 1    angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));    xy = xyrotate(allx,ally,angle);    else        xy = cat(2,allx,ally);        angle = 0;    endexy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);crit = mean(exy(:,1));endif DATA.interactive < 0    [dip, mydip] = GMDip(xy, DATA.energy(1,:),'crit',crit,'label',DATA.idstr);    Cut.gmdprime = mydip.gmdprime;    return;endif ~isfield(E,'sign')    E.sign = 0;endif E.shape(1) == 0    clid = find(DATA.clst(DATA.uid) == DATA.currentcluster);    nid = find(DATA.clst(DATA.uid) ~= DATA.currentcluster);elseif E.sign > 0clid = find(xy(:,1) > crit);nid = find(xy(:,1) <= crit);elseclid = find(xy(:,1) < crit);nid = find(xy(:,1) >= crit);endif replotptshold off;plot(allx,ally,'.','markersize',1);hold on;plot(E.pos([1 3]),E.pos([2 4]),'r-');plot(mean(allx(clid)),mean(ally(clid)),'r+');plot(mean(allx(nid)),mean(ally(nid)),'r+');endcfig = AllV.SetFigure(DATA.tag.hist, DATA,DATA.watcharg{:});subplot(2,1,2);xid = [];showgmfit = 0;fprintf('Cl%d: ',DATA.currentcluster);if isfield(E,'bestcl') && length(E.bestcl) == length(xy) && DATA.usegmcid    fprintf('XY using GM clustering\n');    nid = find(E.bestcl == 1);    clid = find(E.bestcl == 2);    xid=find(E.bestcl ==3);    diffid = find(E.bestcl ==1 & xy(:,1) .* E.sign > crit .* E.sign);    diffcid = find(E.bestcl ==2 & xy(:,1) .* E.sign < crit .* E.sign);elseif E.shape == 0    fprintf('XY using Ellipse\n');    diffcid = [];elseif isfield(E,'gmfit2dman') && DATA.usegmcid    showgmfit = 1;    fprintf('XY using GM 2D clustering\n');    if E.shape == 0        E.bestcl = cluster(E.gmfit2dman,DATA.xy{1});    elseif E.shape == 1        E.bestcl = cluster(E.gmfit2dman,xy);    else        E.bestcl = cluster(E.gmfit2dman,cat(2,allx,ally));    end    nid = find(E.bestcl == 1);    clid = find(E.bestcl == 2);    xid=find(E.bestcl ==3);    diffid = find(E.bestcl ==1 & xy(:,1) .* E.sign > crit .* E.sign);    diffcid = find(E.bestcl ==2 & xy(:,1) .* E.sign < crit .* E.sign);else    diffcid = [];endhold off;plot(xy(nid,1),xy(nid,2),'.','markersize',1);hold on;plot(xy(clid,1),xy(clid,2),'r.','markersize',1);plot(xy(xid,1),xy(xid,2),'g.','markersize',1);if isfield(DATA,'oldclusterpts')    plot(xy(DATA.oldclusterpts{1},1),xy(DATA.oldclusterpts{1},2),'+');endif showgmfit    plot(E.gmfit2dman.mu(1,1),E.gmfit2dman.mu(1,2),'c+','linewidth',2);    plot(E.gmfit2dman.mu(2,1),E.gmfit2dman.mu(2,2),'c+','linewidth',2);endif length(diffcid) && DATA.watchplotsAllV.PlotSpikes(DATA,cat(1,diffid, diffcid),'fixy');figure(cfig);endif plotgm    a = AllV.FitGaussMeans(xy,2);   ezcontour(@(x,y)pdf(a.obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));   if isfield(E,'bestspace') && isfield(E,'mahal')       str = sprintf(' Mahal %.2f (1D %.2f, 2D %.2f)',E.bestspace(2),E.mahal(4),E.mahal(1));   elseif isfield(E,'bestspace')       str = sprintf(' Mahal %.2f (%.2f)',a.mahal,E.bestspace(1));   elseif isfield(E,'mahal')       str = sprintf(' Mahal 1D %.2f, 2D %.2f',E.mahal(4),E.mahal(1));   else       str = sprintf(' Mahal %.2f (%.2f)',a.mahal);   endelseif isfield(E,'bestspace')     str = sprintf('Mahal %.2f 1D %.2f 2D %.2f',E.bestspace(1),E.mahal(4),E.mahal(1));elseif isfield(E,'mahal')    str = sprintf('Mahal %.2f.%.2f',E.mahal(1),E.mahal(4));else    str = [];endif E.shape == 2 || E.space(1) == 6    title(sprintf('ND: %s%s',DATA.gmtypelabels{E.space(2)},str));elseif isempty(pcplot)    title(sprintf('Var-E %s',str));else    title(sprintf('%d: %dvs%d%s',DATA.plottype,pcplot(1),pcplot(2),str));endhdat = get(gcf,'UserData');hdat.elmousept.pos(1) = crit;hdat.elmousept.pos(3) = crit;hdat.elmousept.pos(2) = exy(1,2);hdat.elmousept.pos(4) = exy(2,2);hdat.elmousept.shape = E.shape;if isfield(E,'angle')hdat.elmousept.angle = E.angle;elsehdat.elmousept.angle = 0;endif isfield(E,'cluster')hdat.elmousept.color = DATA.colors{E.cluster+1};elsehdat.elmousept.color = 'r';endhdat.elmousept.down = 0;%hdat.elmousept.h = AllV.DrawEllipse(hdat.elmousept);if E.shape ~= 0     plot(exy(:,1),exy(:,2),'r-');else    hdat.elmousept.h = AllV.DrawEllipse(hdat.elmousept);enddp = (mean(xy(clid,1))-mean(xy(nid,1)))./sqrt(mean([var(xy(clid,1)) var(xy(nid,1))]));hold off;subplot(2,1,1);hold off; if E.shape == 0 && isfield(E,'r')[a,x] = hist(AllV.Rprime(E.r),500);bar(x,a,1);axis('tight');hold on;plot([1 1],get(gca,'ylim'),'r-');else[a,x] = hist(xy(:,1),500);bar(x,a,1);axis('tight');hold on;plot(exy(:,1),get(gca,'ylim'),'r-');endarea = trapz(x,a);if ~isfield(E,'quick')    E.quick = 0;end%[dip, mydip] = FindDip(xy(:,1),DATA.energy(1,:),'eval',crit,'plot','gmix');if plotgm == 2    [dip, mydip] = GMDip(xy, DATA.energy(1,:),'plot','crit',crit,'label',DATA.idstr);    if mydip.converged(1) == 0        AllV.oldFindDip(xy(:,1),DATA.energy(1,:),'eval',crit,'plot',DATA.tag.dips,'gmix');    endelseif isfield(E,'gmfit1d') & strmatch('mu',fieldnames(E.gmfit1d),'exact') & quickmode.quick == 0    if E.shape == 0 && isfield(E,'r')        [dip, mydip] = GMDip(AllV.Rprime(E.r), E.gmfit1d,'crit',1,'label',DATA.idstr);    elseif E.quick == 0  || DATA.quickcutmode.fit1cut; %if its a new cut and gmfit1d is old, don't show        [dip, mydip] = GMDip(xy, E.gmfit1d,'crit',crit,'label',DATA.idstr);    else        dip(1:4) = NaN;        mydip.type = 0;        mydip.dipsize = 0;        mydip.gmdprime = 0;        mydip.cdipsize = 0;    end    mydip.sign = E.sign;elseif quickmode.quick == 0     [dip, mydip] = GMDip(xy, DATA.energy(1,:),'crit',crit,'label',DATA.idstr);else    if E.shape == 1 %manual line, at least cal GM fit starting with crit        [dip, mydip] = GMDip(xy, DATA.energy(1,:),'critonly',crit,'label',DATA.idstr);    end    dip(1:4) = 0;    mydip.sign = 0;    mydip.type = 0;    mydip.dipsize = 0;    mydip.gmdprime = 0;    mydip.cdipsize = 0;endif E.sign == 0    E.sign  = mydip.sign;endmycrit = dip;if isfield(mydip,'gxy') %pdf of GM fit in 1D    plot(mydip.gxy(:,1),sum(mydip.gxy(:,[2 3]),2).*area,'r');    plot(mydip.gxy(:,1),mydip.gxy(:,2).*area,'g');    plot(mydip.gxy(:,1),mydip.gxy(:,3).*area,'g');endif isfield(C,'fitdpparams')    dpxy(:,1) = FitGauss(x,C.fitdpparams(1,:),'eval');    dpxy(:,2) = FitGauss(x,C.fitdpparams(2,:),'eval');    if length(C.fitdprime) > 3        pscale = mean(diff(x))./C.fitdprime(4);    else    pscale = area./trapz(x,sum(dpxy,2));    end        plot(x,dpxy(:,1).*pscale,'m');    plot(x,dpxy(:,2).*pscale,'m');endif mydip.type == 1    plot([dip(1) dip(1)],get(gca,'ylim'),'g-');    plot([dip(2) dip(2)],get(gca,'ylim'),'m-');    plot([dip(3) dip(3)],get(gca,'ylim'),'m--');    if length(dip) > 3    plot([dip(4) dip(4)],get(gca,'ylim'),'c-');    endelseplot([dip(1) dip(1)],get(gca,'ylim'),'g-');plot([dip(2) dip(2)],get(gca,'ylim'),'g-');endif length(dip) >4 && ~isnan(dip(5))plot([dip(5) dip(5)],get(gca,'ylim'),'m-');plot([dip(6) dip(6)],get(gca,'ylim'),'m--');end    if checkdprimes[a,b] = AllV.MaxDprime(xy(:,1));yl = get(gca,'ylim');scale = yl(2)./10;plot(b.crit,abs(b.dps).*scale,'r');endp = AllV.ProbeNumber(DATA);if quickmode.quick == 0    dip = HartigansDipTest(sort(xy(:,1))).*100;    bii = AllV.BimodalCoeff(xy(:,1),1.5);    t = sprintf('P%d%s Dip %.1f(%.1f,%.2f gm%.2f)',p,DATA.probelabel,dip,mydip.dipsize(1),bii,mydip.gmdprime);else    t = sprintf('P%d%s Dip (%.1f gm%.2f)',p,DATA.probelabel,mydip.dipsize(1),mydip.gmdprime);endif isfield(E,'fitdprime')    t = [t sprintf(' G2 %.2f',E.fitdprime(1))];endif isfield(E,'isolation')    t = [t sprintf(' Is %.2f',E.isolation(1))];    if isfield(E,'bestisolation') && isfield(E.bestisolation,'isolation')        t = [t sprintf('(%.2f)',E.bestisolation.isolation(1))];    endendtitle(t);DATA.clid = clid;DATA.nid = nid;Cut = E;Cut.area = area;Cut.dip = mycrit;Cut.angle = angle;Cut.crit = [crit mycrit];Cut.hdip = dip;Cut.gmdprime = mydip.gmdprime;Cut.mydip = [mydip.cdipsize mydip.dipsize];Cut.space = [DATA.plottype pcplot];Cut.shape = E.shape;Cut.y = exy(:,2);DATA.cluster.shape = 1;hdat.cluster = Cut;C = AllV.ClusterFromBoundary(E, Cut);set(gcf,'UserData',hdat);%DATA = AllV.ReplotPCs(DATA,E);if C.shape == 1newE = [];%newE = AllV.BoundaryFromCluster(newE, C, DATA.currentcluster);%DrawLine(newE);end