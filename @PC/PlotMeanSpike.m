function h = PlotMeanSpike(C, p, cluster, varargin)% h = PlotMeanSpike(C, p, cluster, varargin)% or h = PlotMeanSpike(DATA,[ eid probe]);addstr = [];listtype = 'default';linetype = 'default';linew = 2;    plots = [ 1 1];    DATA = [];    h = [];    colors = mycolors('spkcolors');strargs = cell2cellstr(varargin);    j = 1;     while j <= length(varargin)        if isstruct(varargin{j}) && isfield(varargin{j},'ArrayConfig')            DATA = varargin{j};        elseif strncmpi(varargin{j},'addtitile',5)            j = j+1;            addstr = varargin{j};        elseif strncmpi(varargin{j},'autolist',5)            listtype = varargin{j};        elseif strncmpi(varargin{j},'lineonly',5)            plots = [0 1];            linetype = 'chspk';        elseif strncmpi(varargin{j},'imageonly',7)            plots = [1 0];        elseif strncmpi(varargin{j},'oneimage',7) %just image for base cl            plots = [1 0];            cluster = C.cluster;        elseif strncmpi(varargin{j},'smallimage',7)            plots = [4 0];        elseif strncmpi(varargin{j},'exptim',6)            plots = [2 0];        elseif strncmpi(varargin{j},'meanlines',6)            plots = [3 0];        elseif strncmpi(varargin{j},'twoim',5)            plots = [1 2];        elseif strncmpi(varargin{j},'allcluster',10)            plots = [5 0];        end        j = j+1;    end%Cluster structs can have progame too        if isfield(C,'toplevel') && ~isfield(C,'space')         DATA = C;        Clusters = getappdata(DATA.toplevel,'Clusters');        C = Clusters{p(1)}{p(2)};        if length(p) > 2            cluster = p(3);        else            cluster = 1;        end        p = p(2);    end        if isfield(DATA,'probes')        nprobes = DATA.nprobes;    else        nprobes = size(C.MeanSpike.ms,1);    end    if size(C.MeanSpike.ms,1) == 1        chspk = 1;        voff = 0;    else        if isfield(C,'chspk') && size(C.MeanSpike.ms,1) >= max(C.chspk)            chspk = C.chspk;        else            chspk = C.probe(1) + [-1:1];            chspk = chspk (chspk > 0 & chspk <= nprobes);        end        if isfield(DATA,'voffset') && length(DATA.voffset) > max(chspk)            voff = DATA.voffset(chspk)-DATA.voffset(p);        else        voff = [-1:1] .*2;        end    end        if cluster == 0 %plot all        nclusters = 1;        for j = 1:length(C.next)            if isfield(C.next{j},'MeanSpike')                nclusters = nclusters+1;            end        end        if size(C.MeanSpike.ms,1) == 1            PC.PlotMeanSpike(C,1,1,'meanlines',DATA, listtype);            title(sprintf('P%s/%d Ex %.1f Gm %.2f (%.2f) %s',PC.ProbeLabel(p, DATA),cluster,C.exptno,C.mahal(1),C.mahal(2),addstr));            return;        elseif plots(1) == 3            PC.PlotMeanSpike(C,p,1,'meanlines',DATA, listtype);        return;        end        nc = 1;        if sum(plots > 0) > 1            subplot(nclusters,2,1);            np = 2;        else            subplot(nclusters,1,1);            np = 1;        end        if plots(1) == 1            PC.PlotMeanSpike(C, p, 1, 'imageonly',DATA,listtype);        end        if np == 2        subplot(nclusters,2,2);        end        if plots(2) == 2            PC.PlotMeanSpike(C, p, -1, 'imageonly',DATA,listtype);        elseif plots(2) > 0            PC.PlotMeanSpike(C, p, 1, 'lineonly',DATA,listtype);        end        nc = 1;        for j = 1:length(C.next)            if isfield(C.next{j},'MeanSpike')                subplot(nclusters,np,1+nc*np);                if plots(1) == 1                PC.PlotMeanSpike(C, p, j+1, 'imageonly',DATA, listtype);                end                if np == 2                subplot(nclusters,np,2+nc*np);                end                if plots(2)                PC.PlotMeanSpike(C, p, j+1, 'lineonly',DATA, listtype);                end                nc = nc+1;            end        end        return;    end    if sum(plots > 0) > 1        subplot(1,2,1);    end    if cluster > 1        if length(C.next) > cluster-2            C.next{cluster-1}.exptno = C.exptno;            C.next{cluster-1}.exptid = C.exptid;            C = C.next{cluster-1};        end    end    if ~isfield(C,'MeanSpike') %can happen with new cluster        return;    end        if size(C.MeanSpike.ms,1) == 1            p = 1;        elseif p <= 0 && isfield(C,'probe');            p = C.probe(1);        end    if plots(1) == 1 || plots(1) == 4        hold off;        if cluster < 0            h(1) = imagesc(C.MeanSpike.mu);        elseif plots(1) == 4 && isfield(C,'chspk')            h(1) = imagesc(C.MeanSpike.ms(C.chspk,:));        else            h(1) = imagesc(C.MeanSpike.ms);        end        colormap('jet');        if size(C.MeanSpike.ms,1) == 1            p = 1;        elseif p <= 0 && isfield(C,'probe');            p = C.probe(1);        end        line([0 5],[p p],'color','r');        cellstr = '';        cellno = 0;        if ~isempty(DATA)            [~, cellno] = PC.isacell(DATA, C.exptid, p, C.cluster, listtype);            ds = DATA.mahaltype;        else            ds = 'defaultisolation';        end        if abs(cellno) > 0            cellstr = sprintf('Cell%d',cellno);         end        if isfield(C,'dropi')            dropi = C.dropi(3);        else            dropi = [];        end        h = title(sprintf('P%d/%d Ex %.1f I%.1fD%.1f %s %s',p,cluster,C.exptno,PC.DistanceMeasure(C,cluster,ds),dropi,cellstr,addstr));        if ~isfield(DATA,'currentcluster')            set(h,'color', colors{C.cluster+1});                    elseif cluster == DATA.currentcluster            set(h,'color', DATA.colors{cluster+1});        end        set(gca,'xtick',[], 'ytick',[]);        if sum(plots > 0) > 1            subplot(1,2,2);        end    elseif plots(1) == 3 %mean lines        if cluster == 0            subplot(1,1,1);        hold off;        voff(1) = 0;        for j = 1:length(chspk)            mm(j,:,1) = minmax(C.MeanSpike.ms(chspk(j),:));            if cluster == 0            for k = 1:length(C.next)                mm(j,:,k+1) = minmax(C.next{k}.MeanSpike.ms(chspk(j),:));            end            end            mx(j,1) = min(mm(j,1,:));            mx(j,2) = max(mm(j,2,:));            if j > 1                voff(j) = mx(j-1,2)-mx(j,1);            end        end        voff = cumsum(voff).*0.8;                for j = 1:length(chspk)             plot(C.MeanSpike.ms(chspk(j),:)+voff(j),'r','linewidth',2);             hold on;             plot(C.MeanSpike.mu(chspk(j),:)+voff(j),'color',colors{1},'linewidth',2);             if cluster == 0             for k = 1:length(C.next)                 if isfield(C.next{k},'MeanSpike')                     plot(C.next{k}.MeanSpike.ms(chspk(j),:)+voff(j),'color',colors{k+2},'linewidth',2);                 end             end             end        end        else            voff = zeros(size(chspk));            for j = 1:length(chspk)                plot(C.MeanSpike.ms(chspk(j),:)+voff(j),'color',colors{j+1},'linewidth',2);                hold on;            end            if length(chspk) == 1                plot(C.MeanSpike.mu(chspk,:)+voff(j),'color',colors{1},'linewidth',2);            end        end            elseif plots(1) == 5 % plot mean spikes of clusters in all Expts on a probe        if cluster == 1            line(1:length(C.MeanSpike.mu(p,:)),C.MeanSpike.mu(p,:),'Color',...                colors{1},'LineWidth',2)        end        line(1:length(C.MeanSpike.ms(p,:)),C.MeanSpike.ms(p,:),'Color',...            colors{cluster+1},'LineWidth',2)    end    if plots(2)        hold off;        v = std(C.MeanSpike.ms');        if strcmp(linetype,'chspk')            id = chspk';        else            id = find(v > max(v)/2);        end                if length(v) >= p && v(p) < 0.1 %low v range so rescale dp;            x = max(abs(minmax(C.MeanSpike.ms(:))));            dpscale = x/3;        else            dpscale = 1;        end        color = mycolors(8);        for j = 1:length(id)            p = id(j);            h(2) = plot(C.MeanSpike.ms(p,:),'color',colors{j+1},'linewidth',linew);            lh(j) = h(2);            labels{j} = sprintf('%d',p);            hold on;            if isfield(C.MeanSpike,'dp') && size(C.MeanSpike.dp,1) >= j                plot(C.MeanSpike.dp(j,:).*dpscale,'g');            end            plot(C.MeanSpike.mu(j,:),'color',[0.5 0.5 0.5]);        end        if sum(strcmp('nolegend',strargs)) ==0            mylegend(lh,labels,'LowerRight');        end    end