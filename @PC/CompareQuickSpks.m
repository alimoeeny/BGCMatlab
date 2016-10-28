function D = CompareQuickSpks(DATA, varargin)
%%PC.CompareQuickSpks(DATA) Shows QuickSpks for two units currently
%%selecting in cellist
%PC.CompareQuickSpks(DATA, [e p cl; e p cl]) compares specified units
%PC.CompareQuickSpks(DATA, X) where X is an xcorrs strucutre
%PC.CompareQuickSpks(DATA, C) where C cell array of clusters

DATA = GetDataFromFig(DATA);
D.fig = [];
C = []; %CellList
e = [];
A = {};
splittime = []; 
tstrs = {'' ''};
setfits = [0 0];
j = 1;
while j <= length(varargin)
    if isnumeric(varargin{j})
        x = varargin{j};
        e = x(:,1);
        ps = x(:,2);
        clid = x(:,3);
    elseif isfield(varargin{j},'probe')
        e(1:2) = DATA.currentpoint(1);
        ps = floor(varargin{j}.probe);
        clid = round(rem(varargin{j}.probe,1).*10);
    elseif strcmp(varargin{j},'celllist')
        j = j+1;
        C = varargin{j};
    elseif strcmp(varargin{j},'clusters')
        j = j+1;
        clid = varargin{j};
    elseif strcmp(varargin{j},'fits')
        j = j+1;
        setfits = varargin{j};
    elseif strcmp(varargin{j},'splittime')
        j = j+1;
        splittime = varargin{j};
    elseif strcmp(varargin{j},'title')
        j = j+1;
        if ischar(varargin{j})
            tstrs{1} = varargin{j};
        elseif iscell(varargin{j})
            tstrs = varargin{j};
        end
    elseif iscluster(varargin{j})
        A = varargin{j};
        e = CellToMat(A,'exptid');
        ps = CellToMat(A,'probe');
        clid = CellToMat(A,'cluster');
    end
    j = j+1;
end
if isempty(e)
    [icells, e, ps, clid] = PC.SelectedCells(DATA);
end
    
    cells(1) = e(1);
cells(2) = ps(1) +clid(1)./10;
if length(ps) == 1
    cells(3) = cells(2);
else
    cells(3) = ps(2) +clid(2)./10;
end
if isempty(C)
    C = PC.GetValue(DATA, 'CellList');
end
Clusters = PC.GetValue(DATA,'Clusters');
tagb = ['Compare' DATA.tag.spikes];


if ~isempty(A) && isfield(A{1},'clst')
    probelist = unique(cat(2,A{1}.chspk,A{2}.chspk));
elseif length(e) > 1
%Get data for all the clusters    
    A = PC.GetClusterInfo(Clusters,[e(1) ps(1) 1; e(2) ps(2) 1],'clst');
    probelist = unique(cat(2,A{1}.chspk,A{2}.chspk));
else
    A{1} = PC.GetClusterInfo(Clusters,[e(1) ps(1) clid(1)],'clst');
    probelist = A{1}.chspk;
    A{2} = A{1};
    if splittime
        ti = find(A{1}.t > splittime,1);
        istep = max([1 round(ti/500)]);
        spki{1} = 1:istep:ti;
        istep = max([1 round((length(A{1}.t)-ti)/500)]);
        spki{2} = ti:istep:length(A{1}.t);
        F(1) = PC.SetFigure(DATA, DATA.tag.spikes);
        [~, details, Spks{1}] = PC.QuickSpikes(DATA, [e(1) cells(2)],'ispk',spki{1});
        set(get(gca,'title'),'color',DATA.colors{1+clid(1)})
        F(2) = GetFigure(tagb,'parent',DATA.toplevel);
        [~, details, Spks{1}] = PC.QuickSpikes(DATA, [e(1) cells(2)],'ispk',spki{2});
        set(get(gca,'title'),'color',DATA.colors{1+clid(1)})
        PC.SetFigure(DATA, DATA.tag.xyplot);
        PC.PlotClusterXY(DATA, A{1},'ispk',[1:ti],'cellid',clid(1));
        set(get(gca,'title'),'color',DATA.colors{1+clid(1)})

        PC.SetFigure(DATA, DATA.tag.autoxyplot);
        PC.PlotClusterXY(DATA, A{1},'ispk',ti:length(A{1}.t),'cellid',clid(1));
        set(get(gca,'title'),'color',DATA.colors{1+clid(1)})
        return;
    end
end
if length(probelist) ==3 && sum(probelist == [21 23 24]) == 3
    probelist = 21:24;
end
alist = getappdata(DATA.toplevel,'AutoCellList');
qargs = {'showmean'}; %args for quickpsks
ue = unique(e);
if length(ue) ==1  && length(e) > 1 %two cells, one expt, can check simultenaity
    [dups, dupd] = clust.isdup(Clusters{ue}(ps));
    PC.SetFigure(DATA,DATA.tag.xcorr);
    id = find(dups > 0);
    if isempty(id)
        id = find(dups >= 0);
    end
    hold off;
    labels = {};
    for j = id(:)';
        plot(dupd(j).xc);
        [a,b] = ind2sub(size(dups),j);
        labels{end+1} = sprintf('%d<->%d %.2f',a,b,min(dupd(j).e(5:6)));
        hold on;
    end
    if ~isempty(labels)
        legend(labels);
    end   
end
if setfits(1)
    A{1} = clust.SetFit(A{1},setfits(1));
end
F(1) = PC.SetFigure(DATA, DATA.tag.spikes);
[~, details, Spks{1}] = PC.QuickSpikes(DATA, [e(1) cells(2)],'withprobes',probelist,A{1},qargs{:});
DATA.voffset = details.voffset;
X.e = e(1);
X.p = ps(1);
setappdata(F(1),'SpikeData',X);
cella= celllist.get(C,[e(1) ps(1) clid(1)]);
D.title(1) = title(sprintf('Cell%d E%dP%.1f %d spikes%s',cella,e(1),cells(2),A{1}.ncut,tstrs{1}),'color',DATA.colors{clid(1)+1});
PC.SetFigure(DATA, DATA.tag.xyplot);
setappdata(gcf,'SpikeData',X);
PC.PlotClusterXY(DATA, A{1},'cluster',clid(1));
PC.AddCellLabels(DATA, A{1}.exptid, A{1}.probe,'NW');


cellb= celllist.get(C,[e(2) ps(2) clid(2)]);
[F(2), isnew] = GetFigure(tagb,'parent',DATA.toplevel);
if isnew
    hm = uimenu(F(2),'Label','&Spool','Tag','Spool');
    uimenu(hm,'Label','Spool This Probe','Callback',{@PC.PlotMenu, 'spikes', 'Spool'},'tag','spoolspikes');
    uimenu(hm,'Label','Spool Spks and XY','Callback',{@PC.PlotMenu, 'spikes', 'SpoolXY'},'tag','spoolspikesxy');
    uimenu(hm,'Label','Mark','Callback',{@PC.PlotMenu, 'spikes', 'MarkComparison'},'tag','markcomparison');
    bp = [0 0 0.15 0.05];
    uicontrol(F(2),'style','pushbutton','string','->AllVPcs','Units','Normalized','Position',bp,'Tag','FullVButton2',...
            'value',0,'callback',@CallAllVPcs);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(F(2),'style','pushbutton','string','Mark','Units','Normalized','Position',bp,'Tag','MarkButton',...
            'value',0,'callback',{@PC.PlotMenu, 'spikes','MarkComparison'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(F(2),'style','pushbutton','string','NextFit','Units','Normalized','Position',bp,'Tag','NextFitButton',...
            'value',0,'callback',{@PC.PlotMenu, 'spikes','NextFit'});
        
end

if setfits(2)
    A{2} = clust.SetFit(A{2},setfits(2));
end

[~, ~, Spks{2}] = PC.QuickSpikes(DATA, [e(2) cells(3)],'withprobes',probelist,'voffset',details.voffset,'ylim',details.ylim,A{2},qargs{:});

D.fig = F;
X.e = e(2);
X.p = ps(2);
setappdata(F(2),'SpikeData',X);
xcell = NaN;
if cella <0 && cella < 0
    [a,b] = find(squeeze(C(e,:,:)) == -cella);
    truecell = a+b./10;
else 
    truecell = 0;
end

nspk(1) = length(A{1}.times);
nspk(2) = length(A{2}.times);
%Find ohter spikes in A2 that match A1
b = [];
if length(ue) ==1
    [a,~,b] = MatchTimes(A{1}.times,A{2}.t,0.0002);
    [b,c] = Counts(A{2}.clst(a));
    ccounts{2} = b;
    fprintf('Matches for %s onP%d:',str(A{1}),ps(2));
    
    for j = 2:length(b)
        cellx = C(e,ps(2),c(j)-1);
        fprintf(' cl%d(%d) %d',c(j)-1,cellx,b(j));
        if b(j) > nspk(1)./10
            fprintf('*');
            if c(j)-1 ~= clid(2)
                xcl = c(j)-1;
                if alist.I(ue,ps(1),xcl) > 2.7 && alist.D(ue,ps(1),xcl) > 1.5
                    xcell = C(e,ps(1),xcl);
                end
            end
        end
    end
    fprintf('\n');
    
    [b,id] = sort(b,'descend');
    c=c(id)-1; %clusters in A2 matching A1
    
    issue = 'none';
    state = 'ok';
    xstr = ' ';
    if diff(e) ==0
        for j = 1:length(c)
            xstr = sprintf(' %sC%d:%d/%d',xstr,c(j),a(j),ccounts{2}(j));
        end
    else
        for j = 1:length(ccounts{2})
            xstr = sprintf(' %sC%d:%d',xstr,j,ccounts{2}(j));
        end
    end
else
    xstr = '';
end
D.title(2) = title(sprintf('Cell%d E%dP%.1f %d spks%s ',cellb,e(2),cells(3),A{2}.ncut,xstr),'color',DATA.colors{clid(2)+1});





if length(b) > 2 && b(2) > A{1}.ncut./10 % A{1} has two duplicates 
    issue = 'duplicatea';
    if cella < 0 %already marked as duplicate
        if xcell == 0
            state = 'xcell';
        else
            state = 'isdup';
        end
    elseif cella == 0
        state = 'blank';
    elseif cella > 0
        state = 'doublecell';
    end
end

%Find ohter spikes in A1 that match A@
if length(ue) == 1 %Look into duplicate spikes
[b,~,a] = MatchTimes(A{1}.t,A{2}.times,0.0002);
[b,c] = Counts(A{1}.clst(a));
fprintf('Matches for %s on P%d:',str(A{2}),ps(1));
    
for j = 2:length(b)
    cellx = C(e,ps(1),c(j)-1);
    fprintf(' cl%d(%d) %d',c(j)-1,cellx,b(j));
    if b(j) > nspk(1)./10
        fprintf('*');
        if c(j)-1 ~= clid(1)
            xcl = c(j)-1;
            if alist.I(ue,ps(2),xcl) > 2.7 && alist.D(ue,ps(2),xcl) > 1.5
                xcell = C(e,ps(2),xcl);
            end
        end
    end

end
fprintf('\n');

[b,id] = sort(b,'descend');
c=c(id)-1; %clusters in A2 matching A1
if length(b) >1 && b(2) > A{2}.ncut./10 % A{2} has two duplicates 
    if strcmp(issue,'none')
        issue = 'duplicateb';
        if cellb < 0 %already marked as duplicate
            state = 'isdup';
        elseif cellb == 0
            state = 'blank';
        elseif cellb > 0
            state = 'doublecell';
        end
    end
end
end
PC.SetFigure(DATA, DATA.tag.autoxyplot);
setappdata(gcf,'SpikeData',X);
PC.PlotClusterXY(DATA, A{2},'cluster',clid(2));
PC.AddCellLabels(DATA, A{2}.exptid, A{2}.probe,'NW');

PC.SetFigure(DATA,DATA.tag.autocelllist);
x = getappdata(gcf,'ShowDupBoxes');
for j = 1:length(x)
    if myhandle(x(j))
        delete(x(j));
    end
end

showdiffs = 0;
if showdiffs
[ib, ~,ia] = MatchTimes(A{1}.times,A{2}.times',0.0002);
adiff = setdiff(1:nspk(1),ia);
bdiff = setdiff(1:nspk(2),ib);
dups = clust.isdup(A);
%adiff is elements of A1 that are not in A2
figure(F(2));
a = FindMissing(A, Spks, adiff,ib); %other spikes in A2 that match A1(adiff)
PC.PlotSpikes(DATA,cells([1 3]),a, Spks{2}, A{2},'byspk','onecolor',[0 0 0],'holdon','voffset',details.voffset,'notitle');


figure (F(1));
a = FindMissing(A([2 1]), Spks([2 1]), bdiff,ia);
PC.PlotSpikes(DATA,cells([1 2]),a, Spks{1}, A{1},'byspk','onecolor',[0 0 0],'holdon','voffset',details.voffset,'notitle');
end


function CallAllVPcs(src, event)

PC.GuiMenu(src, event, 'CallAllVPcs');


function ispk = FindMissing(A,Spks, adiff, ib)
%Find spikes in A{2} matching A1(adiff)
%adiff are the spikes in A{1} that do not match the spike in A{2}
%so the events in A{2} should NOT be the cluster of hte spike
ispk = [];
if ~isempty(adiff)
    [a,~,b] = MatchTimes(A{1}.times(adiff),A{2}.t',0.0002);
    xspk = a(A{2}.clst(a) ~= A{2}.cluster +1); %only look for spikes not already  in duplicate
    [a,c] = Counts(A{2}.clst(xspk));

    if length(xspk) > 50
        a = ceil(length(xspk)./50);
        xspk = xspk(1:a:length(xspk));
    end
    [ispk, ~, a] = MatchTimes(A{2}.t(xspk),Spks{2}.times./10000,0.0002);
    [a,b] = MatchTimes(A{1}.times(adiff),Spks{2}.times./10000);
    [c,d] = Counts(A{2}.clst(a));
end

