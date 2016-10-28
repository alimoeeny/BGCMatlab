function ploth = PlotSpikes(DATA, pos, spkid, Spks, C, varargin)%ploth = PC.PlotSpikes(DATA, pos, spkid, Spks, C, varargin)dvdt = 1;ploth = [];fixy = 0;labelpos = 'right';colors {1} = [0.5 0.5 0.5];colors {2} = [1 0 0];colors {3} = [0 1 0];if isfield(DATA,'colors')    colors = DATA.colors;else    colors = mycolors('spkcolors');endscale = 1;yl = [];j = 1;quicktest = 1;showax = 1;showtitle = 1;holdoff = 1;showall = 0;showxV =1;fontsiz = 20;voffset = [];[e, p, cl] = pos2ind(pos);plotorder = [];lastcl = 0;onecolor = [];xprobes = [];while j <= length(varargin)    if strncmpi(varargin{j},'colors',6)        j = j+1;        colors = varargin{j};    elseif strncmpi(varargin{j},'fixy',3)        if length(varargin) > j & isnumeric(varargin{j+1})            j = j+1;            yl = varargin{j};        else            fixy = 1;        end    elseif strncmpi(varargin{j},'order',5)        j = j+1;        plotorder = varargin{j};    elseif strncmpi(varargin{j},'ontop',5)        j = j+1;        lastcl = varargin{j};    elseif strncmpi(varargin{j},'empty',5)        showax = 0;        showtitle = 0;    elseif strncmpi(varargin{j},'holdon',5)        holdoff = 0;    elseif strncmpi(varargin{j},'labelpos',8)        j = j+1;        labelpos = varargin{j};    elseif strncmpi(varargin{j},'notitle',7)        showtitle= 0;    elseif strncmpi(varargin{j},'oneprobeonly',7)        showxV = 0;    elseif strncmpi(varargin{j},'onecolor',7) %force one color for these spikes        j = j+1;        onecolor = varargin{j};            elseif strncmpi(varargin{j},'byspkid',5)        DATA.plotspk.bytrial = 0;    elseif strncmpi(varargin{j},'showall',7)        showall = 1;    elseif strncmp(varargin{j},'voffset',6)        j = j+1;        voffset = varargin{j};    elseif strncmp(varargin{j},'withprobe',6)        j = j+1;        xprobes = varargin{j};    end    j = j+1;endif ~isfield(Spks,'values')    fprintf('Missing Spikes for Expt %d\n',e);    return;endif ~isfield(DATA,'plotspk')    DATA.plotspk.bytrial = 0;    DATA.options.usesavedcodes = 0;endfor j = 1:length(Spks)    if isfield(Spks,'maxv') && strcmp(class(Spks(j).values),'int16')        Spks(j).values = double(Spks(j).values) .* double(Spks(j).maxv)./Spks(j).maxint;    endendif fixy    yl = get(gca,'ylim');endchspk = p;if isfield(C,'chspk')    probelist = C.chspk;else    probelist = [];endif ~isempty(xprobes)    probelist = xprobes;endif ~isfield(DATA,'usegmcid')    DATA.usegmcid = 0;endidstr = '';if DATA.plotspk.bytrial && exist('nt','var')    nt = max([DATA.currenttrial 0]);    T = DATA.Expt.Trials(nt);    if isfield(T,'id')        idstr = sprintf('(id%d)',T.id);    endendif isempty(spkid)    for c= chspk(:)';        if c > 0 & c <= DATA.nprobes        for j = 1:2            h = plot([0 30],[0 0],'color',colors{j});        end        hold on;        end    end    if fixy        set(gca,'ylim',yl);    end    if DATA.plotspk.bytrial && showtitle        nt = DATA.currenttrial;        T = DATA.Expt.Trials(nt);        title(sprintf('Ex%d T%d%s %.2f-%.2f: No Spikes  ed%.2f',e,T.Trial,idstr,T.Start(1),T.End(end),T.ed));    end    return;endif spkid(end) > size(Spks(1).values,1)    fprintf('Spike Events Mismatch E%dP%d\n',e,p);    spkid = spkid(spkid <= size(Spks(1).values,1));endif length(spkid) == 0    return;endif isfield(C,'clst') && length(C.clst) >= max(spkid)       clst = C.clst;    if showall        clst(clst < 1) = 1;    endnc = unique(clst(spkid));nc = nc(nc > 0); %excludes excluded trials/timestwoclusters = 1;else    twoclusters = 0;endif isempty(spkid)    ids{2} = DATA.clid;    ids{1} = DATA.nid;    id = ids{2};   nid = ids{1};elseif ~isempty(onecolor)    ids{1} = spkid;    ids{2} = [];    id = ids{1};elseif DATA.usegmcid && length(DATA.gmcid) >= max(spkid)        nc = unique(DATA.gmcid);        for j = 1:length(nc)            id = find(DATA.gmcid(spkid) == nc(j));            ids{j} = spkid(id);        endelseif DATA.options.usesavedcodes    nc = unique(Spks(1).codes(spkid,1));    for j = 1:length(nc)        ids{nc(j)+1} = spkid(find(Spks(1).codes(spkid,1) == nc(j)));    end    id = ids{end};elseif twoclusters    ids = {[] []}; %min necessary    for j = 1:length(nc)    ids{nc(j)} = spkid(find(C.clst(spkid) == nc(j)));    end    id = ids{2};   nid = ids{1};else    nc = unique(Spks(1).codes(spkid,1)+1);    nc = nc(nc > 0);    for j = 1:length(nc)        ids{nc(j)} = spkid(find(Spks(1).codes(spkid,1) == nc(j)-1));    end    if length(ids) ==1        id = ids{1};    else        id = ids{2};    endendispk = p;for j = 1:length(ids)    V{j} = Spks(1).values(ids{j},:)';endif isfield(Spks(1),'xvalues') && showxV    for k = 1:length(Spks(1).xchans)%some of the extra chans are there so that CSD can be calculated. Don't display these unless asked                if isempty(probelist) || ismember(Spks.xchans(k),probelist)            for j = 1:length(ids)                xV{k,j} = Spks(1).xmaxv .* double(squeeze(Spks(1).xvalues(k,:,ids{j})))./(Spks(1).maxint);            end        end    endelse    xV = {};endnspks = length(id);if ~isempty(voffset)    voff = voffset;elseif isempty(xprobes)    voff = DATA.voffset - DATA.voffset(ispk(1));else    voff = DATA.voffset;end    l = size(V{1},1);    if holdoff        hold off;    else        hold on;    end    x = [1:l NaN];if isfield(Spks(1),'Vrange')    dv = max(Spks(1).Vrange)/3;elseif isfield(Spks(1),'VRange')    dv = max(Spks(1).VRange)/3;else    dv = 0;endif strcmp(labelpos,'left')    lbxpos = 1;    halign = 'left';elseif strcmp(labelpos,'right')    lbxpos = size(V{1},1)-2;    halign = 'right';else    lbxpos = size(V{1},1)-2;    halign = 'right';endif ~isempty(xV) && showxV     for j = 1:size(xV,1)         c = Spks(1).xchans(j);         ns = 0;         for k = 1:size(xV,2)            nV = xV{j,k} + voff(c);            ns = ns+length(nV);            if length(nV) == 0            elseif size(nV,1) > 1                 ploth(k,j+1) = plot([0 30],[0 0],'color',colors{k});                nV(l+1,:) = NaN;                set(ploth(k,j+1),'Ydata',reshape(nV,1,prod(size(nV))),'Xdata',repmat(x,1,size(nV,2)));                hold on;            else                ploth(k,j+1) = plot(1:l,nV,'color',colors{j});                hold on;            end            if ~isempty(onecolor) && length(nV) > 0                set(ploth(k,j+1),'color',onecolor);            end         end        if ns > 0            h = text(lbxpos,voff(c)+dv,sprintf('%d',c));            set(h,'fontweight','bold','horizontalalignment',halign, 'fontsize', fontsiz);        end     endendif length(plotorder) == length(V)    sorder = plotorder(:)';else    sorder = 1:length(V);    if lastcl > 0 && lastcl < length(V) %if it equals V, no need to change        sorder(lastcl) = length(V);        sorder(end) = lastcl;    end endfor c= chspk;    if c > 0 & c <= DATA.nprobes        for j = sorder            nV = V{j} + voff(c);            if length(nV) == 0            elseif size(nV,1) > 1                ploth(j,1) = plot([0 30],[0 0],'color',colors{j});                nV(l+1,:) = NaN;                set(ploth(j),'Ydata',reshape(nV,1,prod(size(nV))),'Xdata',repmat(x,1,size(nV,2)));                hold on;            else                ploth(j,1) = plot(1:l,nV,'color',colors{j});                hold on;            end            if ~isempty(onecolor) && length(nV) > 0                set(ploth(j,1),'color',onecolor);            end        end        h = text(lbxpos,voff(c)+dv,sprintf('%d',c));        set(h,'fontweight','bold','horizontalalignment',halign,'fontsize',fontsiz);    endendif isempty(xV) && length(Spks) > 1 %plot additional probes          spkid = spkid(spkid <= length(C.clst));    spkt = round(Spks(1).times(spkid)./10);    for c = 2:length(Spks)        if iscell(Spks)            Spk = Spks{c};        else            Spk = Spks(c);        end        if isfield(Spk,'VRange') && ~isempty(Spk.VRange)            dv = max(Spk.VRange)/3;        elseif isfield(Spk,'Vrange')  && ~isempty(Spk.Vrange)            dv = max(Spk.Vrange)/3;        else            dv = abs(voff(Spk.probe)-voff(Spk.probe))/5;        end        %Spks.times is in 0.1ms ticks        [ix, ida, idb] = intersect(round(Spk.times./10),spkt);        for j = 1:length(nc)            id = find(C.clst(spkid(idb)) ==nc(j));            V{j} = Spks(c).values(ida(id),:)';            l = size(V{1},1);            x = [1:l NaN];            nV = V{j} + voff(Spks(c).probe);            if size(nV,1) > 1                h = plot([0 30],[0 0],'color',colors{j});                nV(l+1,:) = NaN;                set(h,'Ydata',reshape(nV,1,prod(size(nV))),'Xdata',repmat(x,1,size(nV,2)));                hold on;            end        end        h = text(lbxpos,voff(Spks(c).probe)+dv,sprintf('%d',Spk.probe));        set(h,'fontweight','bold','horizontalalignment',halign,'fontsize',fontsiz);    endendif length(yl) == 2    set(gca,'ylim',yl);endif showax == 0    set(gca,'xtick',[],'ytick',[]);endhold off;if showtitleif DATA.plotspk.bytrial    nt = max([DATA.currenttrial 0]);    T = DATA.Expt.Trials(nt);    if isfield(T,'id')        idstr = sprintf('(id%d)',T.id);    else        idstr = '';    endtitle(sprintf('E%d T%d%s %.2f-%.2f: %d/%d ed%.2f',e,T.Trial, idstr,T.Start(1)./10000,T.End(end)./10000, nspks,length(spkid),DATA.Expt.Trials(nt).ed));elsetitle(sprintf('Spikes %d-%d(%.3f-%.3f): %d/%d',...    spkid(1),spkid(end),Spks(1).times(spkid(1)),Spks(1).times(spkid(end)),nspks,length(spkid)));endend