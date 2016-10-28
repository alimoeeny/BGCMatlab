function OptionMenu(a,b, fcn, varargin)    onoff = {'off' 'on'};[DATA, F] = GetDataFromFig(a);srchandle = a;tag = get(a,'tag');if ~strcmp(tag,'PostMenu')    DATA = SaveCallback(DATA,a);end%if this is a button made to mimic a menu, X will have this stuffX = get(a,'UserData');if isfield(X,'menutype') && X.menutype == 1 && ishandle(X.mid)    a = X.mid;endqf = fields(DATA.quickcutmode);autofields = fields(DATA.auto);fullvf = fields(DATA.fullvswitchmode);if strcmp(fcn,'deletetrigger')    DATA.cluster = rmfields(DATA.cluster,'trigset');    SetData(DATA);elseif strcmp(fcn,'loadfromspikes')    DATA.auto.loadfromspikes =  ~DATA.auto.loadfromspikes;    set(a,'checked',onoff{1+DATA.auto.loadfromspikes});    DATA.loadfromspikes = DATA.auto.loadfromspikes;    set(DATA.toplevel,'UserData',DATA);elseif strcmp(fcn,'playfullv')    playdur = varargin{1};    Vall = AllV.mygetappdata(DATA,'Vall');    t = get(gca,'xlim');    nt = max([DATA.currenttrial 1]);    T = DATA.Expt.Trials(nt);    if playdur == 10        t(2) = t(1)+5;        t(1) = t(1)-5;    elseif playdur == 1        t = [T.Start(1) T.End(end)]./10000 +[-0.5 0.5];    elseif playdur == 2        t = [DATA.Expt.Trials(nt-1).Start(1) DATA.Expt.Trials(nt-1).End(end)]./10000 +[-0.5 0.5];    end    if isempty(Vall)        Vall = AllV.Spk2FullV(DATA,t);    end            id = find(Vall.t > t(1)  & Vall.t < t(2));    V = double(Vall.V(DATA.probe(1),id));    V = V./max(abs(V)); %scale to +- 1 for audioplayer    try        P = audioplayer(V,40000);        playblocking(P);    end            elseif strncmp(fcn,'plotspikes',10)    DATA.quickcutmode.(fcn) = ~DATA.quickcutmode.(fcn);    set(a,'checked',onoff{1+DATA.quickcutmode.(fcn)});    set(DATA.toplevel,'UserData',DATA);elseif strcmp(fcn,'postmenu')    if isfield(DATA.guistate,'lastguihandle')         x = DATA.guistate.lastguihandle;        if ishandle(x) && (~isfield(DATA.guistate,'xmenu') || double(x) ~= double(DATA.guistate.xmenu))            if isfield(DATA.guistate,'xmenu') && ishandle(DATA.guistate.xmenu)                delete(DATA.guistate.xmenu);            end            DATA.guistate.xmenu = uimenu(DATA.toplevel,'Label',get(x,'label'),'callback',get(x,'callback'));            SetData(DATA);            AllV.PopupWindow(DATA.toplevel, get(x,'tag'),get(x,'label'),get(x,'callback'));        end    elseif isappdata(DATA.toplevel,'PopupWindowData');        AllV.PopupWindow(DATA.toplevel, '','',[]);    end   SetData(DATA);    elseif strcmp(fcn,'plottedprobes')    AllV.ProbeSelector(DATA, 'plot');elseif strmatch(fcn,'labelwithid')    DATA.plot.labelwithid = ~DATA.plot.labelwithid;    set(a,'checked', onoff{1+DATA.plot.labelwithid});    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,'comment')    DATA.currentprobe = AllV.ProbeNumber(DATA);    PlotComments(DATA.datadir,'parent',DATA.toplevel);%    str  = sprintf('E%dP%d',DATA.exptno,AllV.ProbeNumber(DATA));%    GetString(DATA.tag.comments,DATA.toplevel,{@AllV.AddComment, []},'label',str);    return;elseif strmatch(fcn,'keepspikes')    set(F,'Tag','OldSpikes');        AllV.SetFigure(DATA.tag.spikes, DATA,'front');elseif strmatch(fcn,'keepmeanspikes')    id = unique(DATA.clst);    for j = 2:length(id)        ms{j} = AllV.PlotMeanSpike(DATA,'recalc','cluster',j-1);    end    figure;    voff = DATA.voffset - DATA.voffset(DATA.probe(1));    for j = 2:length(id)        for k = DATA.chspk        plot(ms{j}.ms(k,:)+voff(k),'color',DATA.colors{j},'linewidth',2);        hold on;        end    end    dv = mean(diff(voff(DATA.chspk)))/10;    for k = DATA.chspk        text(1,voff(k)+dv,sprintf('E%dP%d',DATA.exptno,k),'fontsize',DATA.gui.fontsize(1));    endelseif strmatch(fcn,'spoolspikes')    AllV.SetFigure(DATA.tag.spikes, DATA,'front');    AllV.SpoolSpikes(DATA); %ellipse    figure(DATA.toplevel);elseif strmatch(fcn,'autoopt')    flag = get(a,'tag');    DATA.auto.(flag) = ~DATA.auto.(flag);    set(a,'checked',onoff{1+DATA.auto.(flag)});    DATA.checkclusters = DATA.auto.checkcluster;    if DATA.auto.LargeDots        DATA.ptsz = [6 6];    else        DATA.ptsz = [1 1];    end    set(DATA.toplevel,'UserData',DATA);    if sum(strcmp(flag,{'LargeDots'}));        AllV.ReplotPCs(DATA,[]);    end    AllV.SetGUI(DATA);elseif strmatch(fcn,'converttomanual')    DATA.cluster.xyr = DATA.cluster.ellipse.xy;    DATA.cluster.angle = DATA.cluster.ellipse.xy(5);    fprintf('Optimizing Cluster 1...\n');    [DATA.cluster.xyr, details] = clust.OptimizeEllipse(DATA.cluster,'mahaldprimenearsu','xy',DATA.xy{1});    DATA.cluster.angle = DATA.cluster.xyr(5);    DATA.cluster.manual = 3;    DATA.cluster.shape(1) = 0;    DATA.cluster.isolation = details.isolation;    for j = 1:length(DATA.cluster.next)        if isfield(DATA.cluster.next{j},'ellipse')            fprintf('Optimizing Cluster %d...\n',j+1);            [DATA.cluster.next{j}.xyr, details] = clust.OptimizeEllipse(PC.GetClusterInfo(DATA.cluster,j+1,'clst'),...                'mahaldprimenearsu','guess',DATA.cluster.next{j}.ellipse.xy,'xy',DATA.xy{j+1});            DATA.cluster.next{j}.angle = DATA.cluster.xyr(5);            DATA.cluster.next{j}.manual = 3;            DATA.cluster.next{j}.shape(1) = 0;            DATA.cluster.next{j}.isolation = details.isolation;            DATA.cluster.next{j}.triggerset = 0;        end    end    fprintf('Done');    DATA.plottype = 11;    DATA.usegmcid = 0;    DATA = AllV.ClassifyAll(DATA,1,'recluster');    C = PC.GetClusterInfo(DATA.cluster,DATA.currentcluster);    if isfield(C,'xyr')        DATA.elmousept.xyr = C.xyr;        DATA.elmousept.angle = C.angle;    end    AllV.ReplotPCs(DATA,[]);    SetData(DATA);elseif strmatch(fcn,'checkclusters')    DATA.checkclusters = ~DATA.checkclusters;    set(a,'checked',onoff{1+DATA.checkclusters});    DATA.auto.checkclusters = DATA.checkclusters;    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,'quickopts')    flag = get(a,'tag');    DATA.quicksave.(flag) = ~DATA.quicksave.(flag);    set(a,'checked',onoff{1+DATA.quicksave.(flag)});    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,qf)    if ~strcmp(fcn,'quickest')        DATA.quickcutmode.quickest = 0;    end    if strcmp(fcn,'dropi') %turn trigger hist on, in case was turned off by quikest        DATA.quickcutmode.triggerhist = 1;    end    DATA.quickcutmode.(fcn) = ~DATA.quickcutmode.(fcn);    if DATA.quickcutmode.quickest == 1        for j = 1:length(qf)            DATA.quickcutmode.(qf{j}) = 0;        end        DATA.quickcutmode.quickest = 1;    end            c = get(get(a,'parent'),'children');    tags = get(c,'Tag');    for j = 2:length(c)        k = strmatch(get(c(j),'Tag'),qf);        set(c(j),'checked',onoff{1+DATA.quickcutmode.(qf{k})});    end    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,fullvf)    DATA.fullvswitchmode.(fcn) = ~DATA.fullvswitchmode.(fcn);    AllV.SetMenuChecks(get(a,'parent'),DATA.fullvswitchmode);    set(DATA.toplevel,'UserData',DATA);elseif strcmp(fcn,'findbestspaceall')    a = AllV.FindBestSpace(DATA, DATA.cluster,'allclustersedit');elseif strcmp(fcn,'bestspace12')    a = AllV.FindBestSpace(DATA, DATA.cluster,'clusters',2,3,'plotallxy');    DATA.cluster.bestseparation{1} = a;elseif strcmp(fcn,'findbestspace')    a = AllV.FindBestSpace(DATA, DATA.cluster,'plotallxy');    for j = 1:length(a.xnames)        fprintf('Also Good: %s\n',a.xnames{j});    end    if DATA.currentcluster == 1        DATA.cluster.bestisolation = a;    else        DATA.cluster.next{DATA.currentcluster-1}.bestisolation = a;    end    AllV.PlotBestSpace(DATA);elseif regexp(fcn,'[1-9]probefullv')    DATA.plotspk.nfullvprobes = sscanf(fcn,'%d');    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,{'meansummary'})    AllV.PlotAllProbes(DATA, 'spkmean');elseif strmatch(fcn,{'xysummary'})    DATA = AllV.LoadTrigTimes(DATA,1:DATA.nprobes,'savexy');    AllV.SpikeDraw(DATA,[],  'allxy');    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,{'spksummary'})    AllV.SpikeDraw(DATA,[],  'allquickspks');    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,{'plotsummary'})    DATA = AllV.LoadTrigTimes(DATA,1:DATA.nprobes,'savexy');    AllV.SpikeDraw(DATA,[],  'allquickspks');    AllV.SpikeDraw(DATA,[],  'allxy');    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,{'nextfullv'})        args = AllV.NewArgs(DATA,'expt');    AllV.AllVPcs(DATA.toplevel, 'newexpt', DATA.exptno+1, DATA.probeswitchmode,args{:});elseif strmatch(fcn,{'prevfullv'})    args = AllV.NewArgs(DATA,'expt');    AllV.AllVPcs(DATA.toplevel, 'newexpt', DATA.exptno-1, DATA.probeswitchmode);elseif strmatch(fcn,{'newfullv'})    s = inputdlg({'New Expt #' 'Probe #'},'New Expt to load',1,{num2str(DATA.exptno+1) num2str(DATA.probe(1))});    if isempty(s)        return;    end    eid = str2num(s{1});    pid = str2num(s{2});    args{1} = 'tchan';    args{2} = pid;    args{3} = DATA.probeswitchmode;    xargs = AllV.NewArgs(DATA,'expt');    AllV.AllVPcs(DATA.toplevel, 'newexpt', eid, args{:},xargs{:});elseif strmatch(fcn,{'profiling'})    DATA.profiling = ~DATA.profiling;    set(a,'checked',onoff{1+DATA.profiling});    DATA.profiling = 2 * DATA.profiling; %1 = just my timing    set(DATA.toplevel,'UserData',DATA);elseif strmatch(fcn,{'retrigger'})    args = AllV.RetriggerDialog(a,b,'popup');    if ~isempty(args)        fprintf('Retriggering\n');        AllV.AllVPcs(DATA.toplevel,'tchan',AllV.ProbeNumber(DATA), 'needfullv', args{:});    else        fprintf('Retrigger Cancelled\n');            endelseif strmatch(fcn,{'showfullv' 'includeprepost'})    DATA.plotspk.(fcn) = ~DATA.plotspk.(fcn);    set(a,'checked',onoff{DATA.plotspk.(fcn)+1});    set(DATA.toplevel,'UserData',DATA);endfigure(DATA.toplevel); 