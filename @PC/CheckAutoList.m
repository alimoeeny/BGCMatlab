function cFig = CheckAutoList(DATA, varargin)
%PC.CheckAutoList(DATA|fig) Look for merges/swaps. 
%%PC.CheckAutoList(DATA, 'square', [e p cl]) Prints out history of assignments to a
%location
%The list of alignment fixes found by CheckAligment can be modifield with:
%'fixlist' only incudes fixes that were/would be applied
%'shortlist' good fixes that were not applied (? why? Oftern empty)
%'badlist' fix but distance between "matches" > 16
%'midlist' fix but distance between "matches" > 4 < 16
%'nofix' usefix = 0 in alignresult
cFig = 0;

DATA = GetDataFromFig(DATA);
strargs = cell2cellstr(varargin);
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'square',5)
        FollowSquare(DATA, varargin{:});
        return;
    elseif strncmpi(varargin{j},'checksame',5)
        CheckAlignErrs(DATA, varargin{j});
        return;
    end
    j = j+1;
end
CellList = DATA.autolist.CellList;
X = getappdata(DATA.toplevel,'AutoCellList');
setappdata(DATA.toplevel,'AutoCellListBackup',CellList);

%Build List of possible merges based on locations
%where the cell numbe changes in nearby experiments on the same prone
mergepair = [];
for p = 1:size(CellList,2)
    
    C = squeeze(CellList(:,p,:));
    [counts, cells] = Counts(C);
    counts = counts(cells > 0);
    cells = cells(cells > 0);
    replaced = zeros(1,max(cells));
    if sum(counts > 1) > 1
        for k = 2:length(cells)
            for j = 1:k-1
                    cellj = cells(j);
                    cellk = cells(k);
                    [pa,ca] = find(C == cellj);
                    [pb,cb] = find(C == cellk);
                    dp = bsxfun(@minus,pa,pb');
                    [a,b] = find(abs(dp)<2);                   
                    nsame = sum(dp(:) ==0);
                    sameexpts = intersect(pa,pb);
%If the cells are on nearby expts, but not on the same probe, worth a look                    
                    if ~isempty(a) && length(intersect(pa,pb)) < 2
                        mergepair(end+1,:) = cat(1,cells([j k]), pa(a(1)));
                    elseif ~isempty(a)
                        fprintf('Cells %d,%d close but also overlap E%s\n',cells([j k]),sprintf(' %d',sameexpts));
                    end
            end
        end
    end
end

gaps = [];
for j = 1:length(cells)
    [e,p,cl] = celllist.find(CellList,cells(j));
    id = find(diff(e) > 1 & diff(e) < 4);
    if ~isempty(id)
        gaps(end+1).gap = e(id)+1;
        gaps(end).cell = cells(j);
        gaps(end).p = p;
        gaps(end).cl = cl;
    end
end

nrows = 5;
ncols = 4;
if isfield(X,'ratecheck');
    nrows = nrows+1;
end


xscale = 0.95;
[F, isnew] = GetFigure('FixAutoList','parent',DATA.toplevel,'trackpos');
if isnew
    gui.SetDefaults(F);
    outname = [DATA.autodatadir '/' 'MarkList.mat'];
    if exist(outname)
        load(outname);
        setappdata(F,'Markers',marks);
    end
end

ClearPlot(F);
h = uimenu(F,'label','Mark');
uimenu(h,'label','Probe Shift','callback',{@AddMark,'probeshift'});
uimenu(h,'label','Missed Split','callback',{@AddMark, 'split'});
uimenu(h,'label','Bad Split','callback',{@AddMark, 'badsplit'});
uimenu(h,'label','Missed Collpase','callback',{@AddMark, 'collapse'});
uimenu(h,'label','Need Ellipse','callback',{@AddMark, 'ellipse'});
uimenu(h,'label','Bad Match','callback',{@AddMark, 'badmatch'});
uimenu(h,'label','Other','callback',{@AddMark, 'other'});
uimenu(h,'label','Clear All Marks','callback',{@AddMark, 'clearall'});

uimenu(F,'label','Next','callback',@NextMark);

cFig = F;
for j = 1:size(mergepair,1)
    strs{j} = sprintf('%d-%d Ex%d',mergepair(j,:));
end

a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
    'style','text','string','Changes in Cell Number');
b =  uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
    'style','popupmenu','string',strs,'callback',@MergePair);
gui.place(b,'left',a);
h = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
    'style','pushbutton','string','View Spikes','callback',@MergePair);
gui.place(h,'left',b);
rowui = a;

if ~isempty(gaps)
    strs = {};
    for j = 1:size(gaps,1)
        for k = 1:length(gaps(j).gap)
            strs{end+1} = sprintf('Cell %d E%dP%d',gaps(j).cell,gaps(j).gap(k),gaps(j).p(k));
        end
    end
    a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','text','string','Gaps in cell Seq');
    gui.place(a,'up',rowui);
    b =  uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','popupmenu','string',strs,'callback',@MergePair);
    gui.place(b,'left',a);
    h = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','pushbutton','string','View Spikes','callback',@CheckGaps);
    gui.place(h,'left',b);
    rowui = a;
end

if isfield(X.dup,'initialpartiallist')
    x = X.dup.initialpartiallist;
    strs = {};
    for j = 1:size(x,1)
        strs{end+1} = sprintf('E%dP%.1f %.1f',x(j,:));
    end
a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
    'style','text','string','Split Duplicates');
    gui.place(a,'up',rowui);
    rowui = a;
    b =  uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','popupmenu','string',strs,'callback',@ShowDup,'Tag','DuplicateList');
    gui.place(b,'left',a);
    h = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','pushbutton','string','Problems Only','callback',@SwitchDup);
    gui.place(h,'left',b);
    a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','pushbutton','string','Fix','callback',{@ShowDup,'fix'});
    gui.place(a,'left',h);
end

if isfield(X,'ratecheck');
    A = X.ratecheck; 
    id = find([A.maxjump] > 0);
    [row, col] = ind2sub(size(A.maxjump),id);
    strs = {};
    nj = 0;
    for j = 1:length(id);
        x = A.jumps(id(j));
%ratecheck uses real expt nos. Conver to row id here
        for k = 1:length(x.ex(:))
            x.ex(k) = find(DATA.exptlist == x.ex(k));
        end
        for k = 1:size(x.ex,2)
            nj = nj+1;
            fprintf('Cell %d E%d(%d) %.1f\n',row(j),x.ex(:,k),x.sz);
            if diff(x.ex) == 0 %within block jump
                strs{nj} = sprintf('Cell %d E%d %.1f',A.cellid(col(j)),x.ex(2,k),A.jumps(id(j)).sz(k));
            else
                strs{nj} = sprintf('Cell %d E%d-%d %.1f',A.cellid(col(j)),x.ex(1,k),x.ex(2,k),x.sz(k));
            end
            jumpdata(nj,:) = [col(j) row(j) x.ex(:,k)'];
        end
    end
a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
    'style','text','string','Rate Jumps');
    gui.place(a,'up',rowui);
    rowui = a;
    b =  uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','popupmenu','string',strs,'callback',@ShowRateJump,'UserData',jumpdata,'tag','SequenceList');
    gui.place(b,'left',a);
    h = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','pushbutton','string','Next','callback',{@ShowRateJump,'Next'});
    gui.place(h,'left',b);
end

if isfield(X,'selfdup')
    strs = {};
    jumpdata = X.selfdup;    
    for j = 1:size(X.selfdup,1)
        fprintf('Cell %d E%d(%d) %.1f\n',X.selfdup(j,3), X.selfdup(j,1),X.selfdup(j,2));
        strs{j} = sprintf('Cell %d E%d P%d',X.selfdup(j,3), X.selfdup(j,1),X.selfdup(j,2));
    end
    a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','text','string','Duplicates within Probe');
    gui.place(a,'up',rowui);
    rowui = a;
    b =  uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','popupmenu','string',strs,'callback',@ShowSelfDup,'UserData',jumpdata,'tag','SelfDuplicates');
    gui.place(b,'left',a);
    h = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','pushbutton','string','Next','callback',{@gui.NextListItem,'SelfDuplicates'});
    gui.place(h,'left',b);
    a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','pushbutton','string','Fix','callback',{@ShowSelfDup,'fix'});
    gui.place(a,'left',h);
end

if isfield(X,'alignerrs')
    ox = X.alignerrs;
    strs = {};
    for j = 1:length(ox)
        pos(j) = ox(j).pos(1,2) + ox(j).pos(1,1)./100;
        if sum(ismember(min(ox(j).pos(1:2,1)):max(ox(j).pos(1:2,1)),X.shakelist))
            ox(j).shake = 1;
        else
            ox(j).shake = 0;
        end

    end
    [~,id] = sort(pos);    
    x = ox(id);
    shakeexpts = [];
    for j = 1:length(x)
        x(j).orderid = id(j);
        md = max(x(j).matchdistance);
        if x(j).usefix > 0
            if isfield(x,'fixed') && x(j).fixed > 0
                strs{end+1} = sprintf('F:E%dP%.1f %.0f : E%dP%.1f %s  ',x(j).pos(1,1:2),x(j).checkcl,x(j).pos(2,1:2),sprintf('%d,',x(j).matches));
                good(j) = 4;                
            elseif max(x(j).matchdistance) < 4
                strs{end+1} = sprintf('E%dP%.1f %.0f : E%dP%.1f %s  ',x(j).pos(1,1:2),x(j).checkcl,x(j).pos(2,1:2),sprintf('%d,',x(j).matches));
                good(j) = 1;
            else
                strs{end+1} = sprintf('%.1f?E%dP%.1f %.0f : E%dP%.1f %s  ',md,x(j).pos(1,1:2),x(j).checkcl,x(j).pos(2,1:2),sprintf('%d,',x(j).matches));
                good(j) = 2;
                if max(x(j).matchdistance) > 16
                    good(j) = 3;
                end
            end
        else
            if isfield(x,'matches')
                strs{end+1} = sprintf('*E%dP%.0f.%.0f : E%dP%.0f %s  ',x(j).pos(1,1:2),x(j).checkcl,x(j).pos(2,1:2),sprintf('%d,',x(j).matches));
            else
                strs{end+1} = sprintf('*E%dP%.0f.%.0f : E%dP%.0f NoMatch ',x(j).pos(1,1:2),x(j).checkcl,x(j).pos(2,1:2));
            end
            if x(j).aligned == 3
                good(j) = 5;
            else
                good(j) = 0;
            end
        end
        if sum(ismember(min(x(j).pos(1:2,1)):max(x(j).pos(1:2,1)),X.shakelist))
            shake(j) = 1;
        else
            shake(j) = 0;
        end
        ogood(id(j)) = good(j);
    end
    fprintf('\n');
    setappdata(F,'alignlist',id);
    orderlist = id;
    setappdata(F,'goodlist',[good; shake]);
    [a,b] = Counts(good);
    labels = {'nofix', 'short' 'close' 'bad' 'fixed' 'notduplicate' 'badtransition'};
    for j = 1:length(a)
        fprintf('%s:%d\n',labels{b(j)+1},a(j));
    end
    id = [];
    if sum(strcmp('shortlist',strargs))
        id = [id find(good==1)];        
    end
    if sum(strcmp('nofix',strargs))
        id = [id find(good==0)];        
    end
    if sum(strcmp('badlist',strargs))
        id = [id find(good==3)];        
    end
    if sum(strcmp('midlist',strargs))
        id = [id find(good==2)];        
    end
    if sum(strcmp('fixlist',strargs))
        id = [id find(good==4)];        
    end
    if sum(strcmp('notdup',strargs))
        id = [id find(good==5)];        
    end
    id = setdiff(id,find(shake>0));
    if sum(strcmp('shake',strargs))
        id = [id find(shake==1)];        
    end
    if isempty(id)
        id = find(good >= 0 & good < 5 & shake ==0);
    end
    strs = strs(id);
    x = x(id);
    if sum(strcmp('badtransitions',strargs))
        strs = {};
        clear shake;
        x = X.badtransitions;
        CheckTMatches(x);
        for j = 1:length(x)
            strs{j} = sprintf('E%dP%d.%d -E%dP%d.%d',x(j).pos(1,1:3),x(j).pos(2,1:3));
            if sum(ismember(min(x(j).pos(1:2,1)):max(x(j).pos(1:2,1)),X.shakelist))
                shake(j) = 1;
            else
                shake(j) = 0;
            end
            x(j).orderid = j;
%             strs{j} = 'test';
        end
        strs = strs(shake==0);
        x = x(shake==0);
        setappdata(F,'alignlist',1:length(x));
    end
    a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','text','string','Expt Boundary Change');
    gui.place(a,'up',rowui);
    rowui = a;
    b =  uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./ncols xscale./nrows],...
        'style','popupmenu','string',strs,'callback',@ShowBoundary,'Tag','BoundaryList');
    set(b,'UserData',x);
    gui.place(b,'left',a);
    h = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./(ncols+3) xscale./nrows],...
        'style','pushbutton','string','>>','callback',{@gui.NextListItem,'BoundaryList'});
    gui.place(h,'left',b);
    a = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./(ncols+2) xscale./nrows],...
        'style','pushbutton','string','Show','callback',{@ShowBoundary,'showfix'});
    gui.place(a,'left',h);
    b = uicontrol(F,'units','normalized','position',[0.05 0.05 xscale./(ncols+2) xscale./nrows],...
        'style','pushbutton','string','Fix','callback',{@ShowBoundary,'fix'});
    gui.place(b,'left',a);
end

function CheckTMatches(x)

cells = unique([x.cell]);
pos = cat(3,x.pos);
e2 = squeeze(pos(2,1,:));
for j = 1:length(cells)
    id = find([x.cell] == cells(j));
    expts = e2(id);
    [a,b] = sort(expts);
    for c = 1:length(id)
    end
end


function ShowSelfDup(src,b, varargin)
strargs = cell2cellstr(varargin);

line = get(src,'value');
strs = get(src,'string');

DATA = GetDataFromFig(src);
CellList = PC.GetValue(DATA, 'CellList');
FX = getappdata(DATA.toplevel,'AutoCellList');
x = FX.selfdup(line,:);
pos = [x(1) x(2) x(3); x(1) x(5) x(6)];
PC.CompareQuickSpks(DATA,pos);
pos = celllist.find(abs(CellList),x(4),'expts',x(1),'matrix');
FollowSquare(DATA,'square',pos);

function ShowRateJump(src,b, varargin)

strargs = cell2cellstr(varargin);

if sum(strcmp('Next',strargs))
    btn = src;
    src = findobj(allchild(gcf),'flat','tag','SequenceList');
    line = get(src,'value');
    if line < length(src.String)
        line = line+1;
        set(src,'value',line);
    else
        return;
    end
end


line = get(src,'value');
strs = get(src,'string');

DATA = GetDataFromFig(src);
CellList = PC.GetValue(DATA, 'CellList');
%Expts = PC.GetValue(DATA, 'Expts');
X = get(src,'UserData');
ps = X(line,:);
FX = getappdata(DATA.toplevel,'AutoCellList');
GetFigPos(DATA.toplevel);
%exptnames is allexpts, ps is row #
exname = unique(DATA.expnames(DATA.exptlist(ps([3 4]))));
PC.SetFigure(DATA,DATA.tag.rateseq);
[ratecheck, Expts] = PC.PlotClusterRates(DATA, 'rateseqone','expt',exname{1},'cell',ps(1),'nomu');
DATA.currentcell = ps(1);
if diff(ps([3 4])) == 0 %rate change event within a single block
    id = find(ratecheck.exid == ps(3));
    hold off;
    itrial = max([1 ratecheck.errs(id).smjumpt - ceil(ratecheck.errs(id).smw)]);
    X = PlotRateSequence(Expts{id},'color','r','bytime','callback',{@PC.HitTrial, ps(1)});
    smc = smooth(X.rates,ratecheck.errs(id).smw);
    plot(X.times,smc);
    x(1) = X.times(itrial);
    plot([x(1) x(1)], get(gca,'ylim'),'k:');
    t = Expts{id}.Trials(itrial).Start(1);
    
    pos = celllist.find(CellList,ps(1),'expts',ps(3),'matrix');
    PC.CompareQuickSpks(DATA,pos,'splittime',t./10000);
else
    pos = celllist.find(CellList,ps(1),'expts',ps([3 4]),'matrix');
    PC.CompareQuickSpks(DATA,pos);
    PC.CompareSpikeShape(DATA,pos,'drifts',2,'bestdelay');
end
DATA.currentpoint = pos(1,1:2);
DATA.currentcluster = pos(1,3);
SetData(DATA);
FollowSquare(DATA,'square',pos);

function ShowBoundary(a,b, varargin)

F = GetFigure(a);

if ~strcmp(get(a,'tag'),'BoundaryList')
    a = findobj(allchild(F),'flat','Tag','BoundaryList');
end
line = get(a,'value');
strs = get(a,'string');
x = get(a,'UserData');
DATA = GetDataFromFig(a);
GetFigPos(DATA.toplevel);
BoundaryCheck(DATA, F,x(line), line, get(a,'tag'),varargin{:});


function CorrTable(C, xc, x)
bcl = [1:size(xc,2)]+1;
acl = [1:size(xc,1)]+1;
fprintf('    %s',str(C{2},'nocluster'));
for j = 1:length(bcl)
    fprintf(' cl%d ',bcl(j)-1);
end
fprintf('\n');
for j = 1:length(acl)
    fprintf('%s cl%d',str(C{1},'nocluster'),acl(j)-1);
    for k = 1:length(bcl)
        fprintf(' %.2f',xc(j,k));
    end
    fprintf('\n');
end

if sum(x.cldistance) > 0
    for j = 1:length(bcl)
        fprintf(' cl%d ',bcl(j)-1);
    end
    fprintf('\n');
    fprintf('%s cl%d',str(C{1},'nocluster'),acl(1)-1);
    for k = 1:length(x.cldistance)
        fprintf(' %.2f',x.cldistance(k));
    end
    fprintf('\n');
end


function TransitionCheck(DATA, F, x, line,varargin)
strarg = cell2cellstr(varargin);

pos = x.pos;
C = PC.GetClusterInfo(DATA, pos,'clst','addfits');
set(gcbf,'name',sprintf('Checking %s / %s',str(C{1}),str(C{2})));
truepos = pos;
if length(x.cid) > 1
    truepos(1,3) = x.cid(1);
    truepos(2,3) = x.cid(2);
end
CellList = PC.GetValue(DATA.toplevel,'CellList');
cC = PC.GetClusterInfo(DATA, truepos);
fitstr = TransitionLabel(x, C, CellList);

D = PC.CompareQuickSpks(DATA,C,'clusters',x.cid);

GetFigure('CompareSpikes','parent',DATA.toplevel);
PC.CompareSpikeShape(cC,'drifts',2,'bestdelay');
drawnow;
[xc, details] = clust.ShapeCorr(C,'proberange',2,'drift',1,DATA.ArrayConfig);
xmatch = x.M(x.cid(1),x.cid(2),:);
if details.probeshift(x.cid(1),x.cid(2)) ~= 0
    fprintf('Cell %d May have drifted in Expt%d/%d\n',x.cell,pos(1:2,1));
%find 
    %see if htere is a better match on the same probe
    pos(2,2) = pos(1,2);
    A = PC.GetClusterInfo(DATA,pos);
    [M, mdetails] = clust.MahalMatrix(A,'shift');
    minM = min(M,[],3);
    [d(1),match(1)] = min(minM(x.cid(1),:));

    pos = x.pos;
    pos(1,2) = pos(2,2);
    A = PC.GetClusterInfo(DATA,pos);
    [M, mdetails] = clust.MahalMatrix(A);
    minM = min(M,[],3);
    [d(2),match(2)] = min(minM(x.cid(2),:));
    if abs(CellList(x.pos(1,1),x.pos(2,2),match(2))) == x.cell
        goodmatch = 2;
    elseif abs(CellList(x.pos(2,1),x.pos(1,2),match(1))) == x.cell
        goodmatch = 1;
    else
        goodmatch = 0;
    end
    if goodmatch > 0
        fitstr{1} = sprintf('Good Match - probe drift');
    else
        fitstr{1} = sprintf('Probe Drift D%.2f-%.2f',min(minM),max(minM));
    end
elseif xmatch(2) < 4.5 && xmatch(3) < 4.5
    goodmatch = 3;
    fitstr{1} = sprintf('Good Match %s',sprintf(' %.2f',xmatch));
elseif xmatch(1) < 4.5 && xmatch(4) < 4.5
    fitstr{1} = sprintf('Good Match %s',sprintf(' %.2f',xmatch));
    goodmatch = 4;
elseif sum(xmatch < 4.5) > 1
    fitstr{1} = sprintf('Good Match %s',sprintf(' %.2f',xmatch));
    goodmatch = 4;    
end
gui.Annotation(D.fig(2),'clearall');
gui.Annotation(D.fig(1),'clearall');
figure(D.fig(1));
text(2,1,sprintf('Cell%d',x.cell),'color',DATA.colors{x.cid(1)+1},'fontsize',18);
figure(D.fig(2));
text(2,1,sprintf('Cell%d',x.cell),'color',DATA.colors{x.cid(2)+1},'fontsize',18);
if ~isempty(fitstr{1}) 
    figure(D.fig(1));
    h = annotation('textbox',[0.1 0 0.9 0.1],'string',fitstr{1},'linestyle','none','FitBoxToText','On','verticalalignment','bottom');
    gui.Annotation(D.fig(1),'add',h);
end
drawnow;

if x.probeshift ==0
    clust.MahalMatrix(C,'plot',x.cid,'delay',details.timeeshift(x.cid(1),x.cid(2)));
    minM = squeeze(min(x.M,[],3));
    maxM = squeeze(max(x.M,[],3));
    [a, ai] = min(minM(x.cid(1),:));
    [b, bi] = min(maxM(x.cid(1),:));
    fitstr{2} = sprintf('Best Mathces are %d,%d, currently %d\n',ai,bi,x.cid(2));
    clust.MahalMatrix(C,x.M,'cl',x.cid(1)); %prints table
else
    fitstr{2} = sprintf('Cell Drifts %d probes');
end

if ~isempty(fitstr{2}) 
    figure(D.fig(2));
    h = annotation('textbox',[0.1 0 0.9 0.1],'string',fitstr{2},'linestyle','none','FitBoxToText','On','verticalalignment','bottom');
    gui.Annotation(D.fig(2),'add',h);
end
set(gcbf,'name',sprintf('FixAutoList %d/%d',line,length(getappdata(gcbf,'alignlist'))));

function BoundaryCheck(DATA, F, x, line,listtag,varargin)
strarg = cell2cellstr(varargin);

if isfield(x,'M')
    TransitionCheck(DATA,F,x, line,listtag, varargin{:});
    return;
end

Expts = getappdata(DATA.toplevel,'Expts');
pos = x.pos;
C = PC.GetClusterInfo(DATA, pos,'clst','addfits');
CorrTable(C, x.xc, x);
CheckData = getappdata(F,'CheckData');

fitstr = '';
CellList = PC.GetValue(DATA, 'CellList');
setappdata(DATA.toplevel,'LastCellList',CellList);

xlist = getappdata(F,'alignlist');
strwin = 0;
clid = [C{1}.cluster C{2}.cluster];
x.clusterfixed = 0;
if isfield(C{1},'fixed') && C{1}.fixed > 0
    x.clusterfixed = 2;
elseif isfield(C{2},'fixed') && C{2}.fixed > 0
    x.clusterfixed = 2;
end
FX = getappdata(DATA.toplevel,'AutoCellList');
if isfield(FX,'fixit')
    if FX.fixit(pos(1,1),pos(1,2)) > 0
        x.clusterfixed = 1;
    elseif FX.fixit(pos(1,1),pos(1,2)) > 0
        x.clusterfixed = 1;
        
    end
end
if isfield(x,'duplicated')
    fprintf('E%dP%d %d,%d duplicates %d in E%dP%d id%d\n',x.pos(2,1:2),x.duplicates(1:2), x.duplicated,x.pos(1,1:2),x.orderid);
%show cluster that has different id. Dangerous. Might want to look at xc
%values
%C{2} is the one that seems to be split.  CHeck its not a drift probelm
    Expt = Expts{pos(2,1)};
    [~, dcheck] = clust.Check(C{2},'drift',FX.ArrayConfig,Expt);
    clid(1) = x.duplicated;
    clid(2) = x.duplicates(1);
    if dcheck.driftswap > 0
        fprintf('Drift %d (%.2f) in %s\n',dcheck.driftswap,max(dcheck.flipscores),str(C{2}));
        x.driftswap = dcheck.driftswap;
    end
   [fitstr, strwin] = FitLabel(x, C, CellList);            
    fprintf('Fix %d(%d). D %.2f,%.2f Aligned:\n',x.usefix,x.fixed,x.matchdistance,x.aligned);
    xp = setdiff(x.duplicates,x.duplicated);
    ok = clust.GoodFix(x, CellList);
    pos(2,3) = xp(1);
end

if ~isfield(CheckData,'showfix')
    CheckData.showfix = 0;
end
if CheckData.showfix == 1 %'Show' Button toggles
    CheckData.showfix = 0;
    strarg = setdiff(strarg,'showfix');
    Cl = getappdata(DATA.toplevel,'LastCellList');
    PC.SetValue(DATA.toplevel,'CellList',Cl);

end

[fiststr, strwin] = FitLabel(x, C, CellList);
if sum(strcmp('showfix',strarg)) && x.usefix > 0
    CheckData.showfix = 1;
    fit = x.newfit(x.usefix);
    if x.reversed
        C = C([2 1]);
    end
    if isfield(x,'driftswap') && x.usefix ~=1 
        cprintf('red','%s looks like spike drift, but best resoluation was not collapse',str(C{2}));
        C{2} = clust.SetFit(C{2},[],'combine',[dcheck.cla dcheck.clb]);
        C{2}.eckercluster.SU(dcheck.clb) = [];
        fitstr{2} = sprintf('Merging %d with %d in %s',dcheck.cla,dcheck.clb,str(C{2}));
        fitstr{1} = sprintf('Drifting waveform makes %d <-> %d',x.duplicates(1:2));
    elseif x.usefix ==1 %merge C2
        C{2} = clust.SetFit(C{2},fit.fitnumber);
        cl = x.duplicates;
        d(1) = C{2}.eDistance(cl(1)+1,cl(2)+1);
        d(2) = C{2}.eDistance(cl(2)+1,cl(1)+1);
        fitstr{1} = sprintf('Alt Fit %d Merging Clusters %s in E%dP%d',x.usefix,sprintf('%d ',x.duplicates),x.pos(2:1:2));
        strwin = 2;
        newCellList = celllist.ApplyFix(CellList, x, x.newfit(x.usefix));
        setappdata(DATA.toplevel,'LastCellList',CellList);
        PC.SetValue(DATA.toplevel,'CellList',newCellList);
    elseif x.usefix == 2 %split C1
        fitstr{1} = sprintf('Using Alternate Fit %d Splitting Cluster %d in E%dP%d. %d and %d',x.usefix,x.duplicated,x.pos(1,1:2),fit.newmatches(1:2));
        strwin = 1;
        C{1} = clust.SetFit(C{1},fit.fitnumber);        
        [newC, dcheck] = clust.Check(C{1},'drift');
        [newCellList, D] = celllist.ApplyFix(CellList, x, x.newfit(x.usefix));
        setappdata(DATA.toplevel,'LastCellList',CellList);
        DATA = PC.SetValue(DATA.toplevel,'CellList',newCellList);
        if isfield(dcheck,'driftswap') && dcheck.driftswap == 1
            fitstr{1} = sprintf('%s !drift%d',fitstr{1},dcheck.driftswap);
        end
        if length(x.splitcl) == length(fit.newmatches)
            for j = 1:length(x.splitcl)
                if fit.newmatches(j) > 0
                oldcell(j,1) = CellList(fit.pos(1),fit.pos(2),fit.newmatches(j));
                oldcell(j,2) = CellList(x.pos(2,1),x.pos(2,2),x.splitcl(j));
                end
            end
        end
        skip = sum(~ismember(oldcell(:,1),[oldcell(:,2); 0]));
        if skip == 0 && ~isfield(D,'celltype') %celllist now set by ApplyFix
            for j = 1:length(x.splitcl)
                CellList(fit.pos(1),fit.pos(2),fit.newmatches(j)) = oldcell(j,2);
            end
            DATA.autolist.CellList = CellList;
            SetData(DATA);

%            PC.SetValue(DATA,'CellList',CellList);
        end
    else
        fitstr{2} = sprintf('Best fit is new ellipse');
%        C{1} = clust.SetFit(C{1},fit.fitnumber);
    end
    fprintf('%s\n%s\n',fitstr{1},fitstr{2});
    GetFigure('CompareSpikes','parent',DATA.toplevel);
    PC.CompareSpikeShape(C,'drifts',2,'bestdelay');

elseif sum(strcmp('fix',strarg))
    [newC, newX] = clust.AlignAutoFits(C,'checkcl',x.duplicated,'noellipse','spaceone');
    if newX.reversed
        C = C([2 1]);
    end
    if newX.aligned == 1 || newX.aligned == -1
        C = C([2 1]);
        [newC, newX] = clust.AlignAutoFits(C,'checkcl',x.duplicated,'noellipse','spaceone');
    end
    [fitstr, strwin] = FitLabel(newX,C, CellList);
        
    fixes = find(newX.dosplit > 0);
    if length(fixes) > 1
        fixstrs{1} = sprintf('Collpase E%d',C{2}.exptno);
        fixstrs{2} = sprintf('Split E%d Fit',C{1}.exptno);
        fixstrs{3} = sprintf('Split E%d Ellipse',C{1}.exptno);
        for j = 1:length(fixes)
            if ~isempty(newX.newfit(fixes(j)).mahal)
            end
        end
        if newX.aligned > 0
            xstr = sprintf(' But aligned');
        elseif newX.usefix == 0
            xstr = sprintf(' But not  used');
        else
            xstr = '';
        end
           
            
        if newX.bestfix > 0
            q = sprintf('Best is %s%s',fixstrs{newX.bestfix},xstr);
            fixstrs = {fixstrs{fixes} fixstrs{newX.bestfix}};
        else
            q = 'Show Which Fix?';
            fixstrs = fixstrs(fixes);
        end
        yn = gui.Dlg(q,DATA.toplevel, fixstrs);
        usefix = find(strcmp(yn,fixstrs),1);
        nofix = setdiff(fixes,fixes(usefix));
        newX.dosplit(nofix) = 0;
        newX.usefix = fixes(usefix);
    end
 
    if newX.aligned > 0
        if newX.usefix ==2 
%Split of C{1} does better matching of clusters in C{2} even though there
%were no true duplications. 
            C{1} = clust.SetFit(C{1},newX.newfit(2));
        else
            fitstr{1} = sprintf('%s. Not Really',fitstr{1});
        end
        
    elseif newX.usefix == 0
        fitstr{1} = sprintf('%s. No Change',fitstr{1});
    elseif newX.dosplit(2) > 0 %split C{1} using existing fit
        C{1} = clust.SetFit(C{1},newX.newfit(2));
    elseif newX.dosplit(3) > 0 %split C{1} with new ellipses
        C{1} = clust.SetFit(newC{1},newX.newfit(3));        
    elseif newX.dosplit(1) > 0 
        C{2} = clust.SetFit(C{2},newX.dosplit(1));
    elseif newX.aligned > 0
        CorrTable(C,newX.xc,newX);
        [newC, newX] = clust.AlignAutoFits(C([2 1]),'checkcl',x.duplicated,'noellipse','spaceone');
    end
    if newX.usefix > 0 && newX.aligned ==0
        [newCellList, details] = celllist.ApplyFix(CellList, newX, newX.newfit(newX.usefix));
        [newCellList, details] = CheckChanges(newCellList,CellList, details, C);
        setappdata(DATA.toplevel,'LastCellList',CellList);
        DATA = PC.SetValue(DATA.toplevel,'CellList',newCellList);
    end
else
%can only do this for existing clusters    
if x.usefix == 2
    elseif x.usefix == 1
    else
    GetFigure('CompareSpikes','parent',DATA.toplevel);
    PC.CompareSpikeShape(DATA,pos,'drifts',2,'bestdelay');
    fprintf('%s Cl%d is duplicated in %s Cl%s\n',str(C{1}),x.duplicated,str(C{2}),sprintf('%d ',x.duplicates));
    end
end
PC.PlotCellList(DATA,'autolist');
GetFigure(DATA.tag.autocelllist);
PC.ClearSelections(DATA,0,0);
PC.MarkSelections(DATA, C);

D = PC.CompareQuickSpks(DATA,C,'clusters',clid);
gui.Annotation(D.fig(2),'clearall');
gui.Annotation(D.fig(1),'clearall');
SetMark(F,x, listtag);

if ~isempty(fitstr{1}) 
    figure(D.fig(1));
    h = annotation('textbox',[0.1 0 0.9 0.1],'string',fitstr{1},'linestyle','none','FitBoxToText','On','verticalalignment','bottom');
    gui.Annotation(D.fig(1),'add',h);
end
if length(fitstr) > 2 && ~isempty(fitstr{3})
    figure(D.fig(1));
    h = text(2,1,fitstr{3},'color',DATA.colors{1+x.duplicated},'FontSize',21);
end
if ~isempty(fitstr{2}) 
    figure(D.fig(2));
    h = annotation('textbox',[0.1 0 0.9 0.1],'string',fitstr{2},'linestyle','none','FitBoxToText','On','verticalalignment','bottom');
    gui.Annotation(D.fig(2),'add',h);
    tx = 5;
    if length(fitstr) > 3 && ~isempty(fitstr{4})
        h = text(2,1,fitstr{4},'color',DATA.colors{1+x.duplicates(1)},'FontSize',21);
        sz = get(h,'extent');
        tx = sz(1)+sz(3);
    end
    if length(fitstr) > 4 && ~isempty(fitstr{5})
        h = text(tx,1,fitstr{5},'color',DATA.colors{1+x.duplicates(2)},'FontSize',21);
    end
end


figure(D.fig(2));
Cl = PC.GetValue(DATA,'CellList');
xp = 0;
e = pos(1,1);
p = pos(1,2);
for j = 1:length(x.duplicates)
    xstr = sprintf('%d ',x.duplicates(j));
    clid = x.duplicates(j);
    if Cl(e,p,clid) > 0
        xstr = sprintf('%sC%d',xstr, Cl(e,p,clid));
    elseif Cl(e,p,clid) < 0
        xstr = sprintf('%sD%d',xstr, -Cl(e,p,clid));
    end
    h = text(xp,1.0,xstr,'units','normalized','FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom');
    set(h,'color',DATA.colors{clid+1});
    te = get(h,'extent');
    xp = xp + te(3);
    if clid ~= C{1}.cluster
        xstr = sprintf('%d',sum(C{1}.clst == clid));
        te = get(D.title(2),'extent');
        text(te(1)+te(3),te(2),xstr,'verticalalignment','bottom','color',DATA.colors{clid+1},'HorizontalAlignment','left','fontweight','bold');
    end
end
figure(D.fig(1));
e = pos(2,1);
p = pos(2,2);
clid = x.duplicated;
xstr = sprintf('%d ',x.duplicated);
    if Cl(e,p,clid) > 0
        xstr = sprintf('%sC%d',xstr, Cl(e,p,clid));
    elseif Cl(e,p,clid) < 0
        xstr = sprintf('%sD%d',xstr, -Cl(e,p,clid));
    end
    h = text(0,1.0,xstr,'units','normalized','FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom');
    set(h,'color',DATA.colors{x.duplicated+1});
    
DATA.currentpoint = pos(1,1:2);
DATA.currentcluster = pos(1,3);
SetData(DATA);
setappdata(F,'CheckData',CheckData);
PC.SetValue(DATA.toplevel,'LastCellList');

function [C, D] = CheckChanges(C, O, D, Clusters)

if isfield(D, 'changes')
    for j = 1:size(D.changes,1)
        x = PC.GetClusterInfo(Clusters, D.changes(j,1:3));
        if x.isolation(1) < 2.8 || x.dropi(3) < 1.8 && D.changes(j,4) > 0
            C(D.changes(j,1),D.changes(j,2),D.changes(j,3)) = 0;
            D.changes(j,4) = -D.changes(j,4);
            D.rejectchange(j,:) = [x.isolation(1) x.dropi(3)];
            cprintf('blue','Undoing Change at E%dP%d.%d Cell%d\n',D.changes(j,:));
        end
    end
end

function mark = SetMark(F, x, tag)

mark.pos = x.pos;
if isfield(x,'usefix')
    mark.type = 'align';
else
    mark.type = '?';
end
mark.tag = tag;
mark.time = now;
mark = CopyFields(mark,x,'orderid','cells','matchcells');
setappdata(F,'CurrentMarker',mark);


function NextMark(src, event, reason)
F = GetFigure(src);
DATA = GetDataFromFig(F);
marks = getappdata(F,'Markers');
markcount = getappdata(F,'markcounter');
if isempty(markcount) || markcount == 0
    markcount = 1;
else
    markcount = markcount+1;
end
if markcount > length(marks)
    markcount = 1;
end
setappdata(F,'markcounter',markcount);
mark = marks{markcount};
if isfield(mark,'tag')
    list = findobj(F,'tag',mark.tag);
    items = get(list,'UserData');
    id = find([items.orderid] == mark.orderid);
    x = items(id);
    set(list,'value',id);
    ShowBoundary(list,[]);
    if isfield(mark,'matchcells')
        fprintf('Cells at list time were %s\n',sprintf('%d ',mark.matchcells));
    end
end
fprintf('Marked %d %d/%d Marks',mark.orderid,markcount,length(marks));
%PC.CompareQuickSpks(DATA, mark.pos);

function AddMark(src, event, reason)

F = GetFigure(src);
DATA = GetDataFromFig(F);
if strcmp(reason,'clearall')
    rmappdata(DATA.toplevel,'Markers');
    return;
end
mark = getappdata(F,'CurrentMarker');
mark.reason = reason;
findcell.AddMark(DATA, mark);

function fitstr = TransitionLabel(x, C, CellList)

if sum(x.clusterok) < 2
    xcell = 0;
    fitstr{1} = sprintf('Combining fit would loose cell%d')
else
    cells(1) = CellList(x.pos(1,1),x.pos(1,2),x.cid(1));
    cells(2) = CellList(x.pos(2,1),x.pos(2,2),x.cid(2));
    fitstr{1} = sprintf('Distances %s',sprintf(' %.1f',x.mahal));
    fitstr{2} = 'supicious transition';
    fitstr{3} = sprintf('%dCell%d',x.cid(1),cells(1));
    fitstr{4} = sprintf('%dCell%d',x.cid(1),cells(1));
    fitstr{5} = sprintf('%dCell%d',x.cid(2),cells(2));
end

function [fitstr, strwin] = FitLabel(x, C, CellList)

if isfield(x,'duplicated')
    cells(1) = CellList(x.pos(1,1),x.pos(1,2),x.duplicated);
    cells(2) = CellList(x.pos(2,1),x.pos(2,2),x.duplicates(1));
    cells(3) = CellList(x.pos(2,1),x.pos(2,2),x.duplicates(2));
elseif isfield(x,'cid')
    cells(1) = CellList(x.pos(1,1),x.pos(1,2),x.cid(1));
    cells(2) = CellList(x.pos(2,1),x.pos(2,2),x.cid(2));
end

if x.aligned == 3
    if isfield(x,'cldistance')
        fitstr{2} = sprintf('%d(%.2f,%.1f),%d(%.2f,%.1f)Not both Duplicates of %d',x.duplicates(1),x.xc(1),x.cldistance(1),x.duplicates(2),x.xc(2),x.cldistance(2),x.duplicated);
    else
        fitstr{2} = sprintf('%d(%.2f),%d(%.2f)Not both Duplicates of %d',x.duplicates(1),x.xc(1),x.duplicates(2),x.xc(2),x.duplicated);
    end    
    if isfield(x,'cell')
        fitstr{3} = sprintf('Cell %d',x.cell);
    else
        fitstr{3} = sprintf('Cell %d',cells(1));
    end
    fitstr{4} = sprintf('%dCell%d',x.duplicates(1),cells(1));
    fitstr{5} = sprintf('%dCell%d',x.duplicates(2),cells(3));
elseif isfield(x,'duplicates')
    fitstr{2} = sprintf('%d,%d duplicates %d.',x.duplicates(1:2),x.duplicated);
    if isfield(x,'xc') && length(x.xc) > 1
        if isfield(x,'cldistance')
            fitstr{2} = sprintf('%d(%.2f,%.1f),%d(%.2f,%.1f) duplicates %d.',x.duplicates(1),x.xc(1),x.cldistance(1),x.duplicates(2),x.xc(2),x.cldistance(2),x.duplicated);
        else
            fitstr{2} = sprintf('%d(%.2f),%d(%.2f) duplicates %d.',x.duplicates(1),x.xc(1),x.duplicates(2),x.xc(2),x.duplicated);
        end
    end
    fitstr{3} = sprintf('%d,%d <-> %d',x.duplicates(1:2),x.duplicated);
elseif x.aligned == 1
    fitstr{2} = 'OK';
elseif x.aligned > 1
    fitstr{2} = 'Similar but OK';
end
if isfield(x,'clusterfixed')
    if x.clusterfixed == 2
        fitstr{2} = [fitstr{2} '*'];
    elseif x.clusterfixed ==1
        fitstr{2} = [fitstr{2} '(*)'];
    end
elseif sum(CellToMat(C,'fixed'))
    fitstr{2} = [fitstr{2} '*'];
end
    
fitstr{1} = '';

if x.usefix == 0
    if x.aligned == 2
        fitstr{2} = sprintf('%s. %d(%.2f) matches. %d(%.2f)',fitstr{2},x.match,x.cldistance(1),x.duplicates(2),x.cldistance(2));
    else
        fitstr{2} = sprintf('%s. No Fix',fitstr{2});
    end
    strwin = 1;
elseif x.usefix == 2
    strwin = 1;
    fitstr{1} = sprintf('Split with Fit %.1f',x.newfit(2).fitnumber);
elseif x.usefix == 1
    strwin = 1;
    fitstr{2} = sprintf('%s Collapse with Fit %.1f',fitstr{2},x.newfit(1).fitnumber);
else
    strwin = 1;
    fitstr{2} = sprintf('%s Usefix%d',fitstr{2},x.usefix);
end
if isfield(x,'driftswap') && x.usefix ~=1
        fitstr{2} = sprintf('Drifting waveform makes %d <-> %d',x.duplicates(1:2));
end

function SwitchDup(src,b)

s = get(src,'String');
if strcmp(s,'Problems Only')
    set(src,'String','All');
else
    set(src,'String','Problems Only');
end



function ShowDup(a,b, varargin)


if ~strcmp(get(a,'tag'),'DuplicateList')
    F = GetFigure(a);
    a = findobj(allchild(F),'flat','Tag','DuplicateList');
end

line = get(a,'value');
strs = get(a,'string');
DATA = GetDataFromFig(a);
GetFigPos(DATA.toplevel);
response = DuplicateCheck(DATA, strs{line},varargin{:});
while ~isempty(response)
    if strcmp(response,'Next')
        line = line+1;
        response = DuplicateCheck(DATA, strs{line});
    else
        response = '';
    end
end


function response = DuplicateCheck(DATA, s, varargin)
strargs= cell2cellstr(varargin);
cells = sscanf(s,'E%dP%f %f');
e = floor(cells(1));
p = floor(cells(2:3));
clid = round(rem(cells(2:3),1).*10);
C = PC.GetValue(DATA, 'CellList');
Clusters = PC.GetValue(DATA,'Clusters','withfits');

A = PC.GetClusterInfo(Clusters,[cells(1) cells(2); cells(1) cells(3)],'clst');
probelist = unique(cat(2,A{1}.chspk,A{2}.chspk));
alist = getappdata(DATA.toplevel,'AutoCellList');

[dups, dupd] = clust.isdup(Clusters{cells(1)}(floor(cells(2:3))));

if sum(strcmp('fix',strargs))
    [newC, newX] = clust.AlignAutoFits(A);
end
PC.SetFigure(DATA,DATA.tag.xcorr);
id = find(dups > 0);
hold off;
labels = {};
for j = id(:)';   
    plot(dupd(j).xc);
    [a,b] = ind2sub(size(dups),j);
    labels{end+1} = sprintf('%d<->%d',a,b);
    hold on;
end
legend(labels);
    


F(1) = PC.SetFigure(DATA, DATA.tag.spikes);
[~, details, Spks{1}] = PC.QuickSpikes(DATA, cells([1 2])','withprobes',probelist);
DATA.voffset = details.voffset;

cella = celllist.get(C,[cells(1) cells(2)]);
cl = round(rem(cells(2:3),1).*10);
title(sprintf('Cell%d E%dP%.1f %d spikes',cella,cells(1),cells(2),A{1}.ncut),'color',DATA.colors{cl(1)+1});
PC.SetFigure(DATA, DATA.tag.xyplot);
PC.PlotClusterXY(DATA, A{1});

tagb = ['Compare' DATA.tag.spikes];
cellb = celllist.get(C,cells([1 3])');
F(2) = GetFigure(tagb,'parent',DATA.toplevel,'trackpos');
[~, ~, Spks{2}] = PC.QuickSpikes(DATA, [cells(1) cells(3)],'withprobes',probelist,'voffset',details.voffset,'ylim',details.ylim);

xcell = NaN;
if cella <0 && cella < 0
    [a,b] = find(squeeze(C(e,:,:)) == -cella);
    truecell = a+b./10;
else 
    truecell = 0;
end

[a,b] = find(squeeze(C(e,:,:)) ==16);
nspk(1) = length(A{1}.times);
nspk(2) = length(A{2}.times);
%Find ohter spikes in A2 that match A1
[a,~,b] = MatchTimes(A{1}.times,A{2}.t,0.0002);
[b,c] = Counts(A{2}.clst(a));
fprintf('Matches for %s onP%d:',str(A{1}),p(2));
    
for j = 2:length(b)
    cellx = C(e,p(2),c(j)-1);
    fprintf(' cl%d(%d) %d',c(j)-1,cellx,b(j));
    if b(j) > nspk(1)./10
        fprintf('*');
        if c(j)-1 ~= clid(2)
            xcl = c(j)-1;
            if alist.I(e,p(1),xcl) > 2.7 && alist.D(e,p(1),xcl) > 1.5
                xcell = C(e,p(1),xcl);
            end
        end
    end
end
fprintf('\n');

[b,id] = sort(b,'descend');
c=c(id)-1; %clusters in A2 matching A1

issue = 'none';
state = 'ok';
if c(2) > 0
    xstr = sprintf('C%d %d',c(2),a(2));
else
    xstr = '';
end

title(sprintf('Cell%d E%dP%.1f %d spks. c%d %d match%s ',cellb,cells(1),cells(3),A{2}.ncut,c(1),a(1),xstr),'color',DATA.colors{cl(2)+1});
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

[b,~,a] = MatchTimes(A{1}.t,A{2}.times,0.0002);
[b,c] = Counts(A{1}.clst(a));
fprintf('Matches for %s on P%d:',str(A{2}),p(1));
    
for j = 2:length(b)
    cellx = C(e,p(1),c(j)-1);
    fprintf(' cl%d(%d) %d',c(j)-1,cellx,b(j));
    if b(j) > nspk(1)./10
        fprintf('*');
        if c(j)-1 ~= clid(1)
            xcl = c(j)-1;
            if alist.I(e,p(2),xcl) > 2.7 && alist.newD(e,p(2),xcl) > 1.5
                xcell = C(e,p(2),xcl);
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

if strcmp(issue,'none') %real dup calc did not match up with the counts here
    dupb = find(dups(:,A{2}.cluster)>0);
    dupa = find(dups(A{1}.cluster,:)>0);
    if length(dupa) > 1
        pid = 1;
        xcl = setdiff(dupa, A{1}.cluster);
        issue = 'duplicatea';
    elseif length(dupb) > 1
        pid = 2
        xcl = setdiff(dupb, A{2}.cluster);
        issue = 'duplicateb';
    else %should' happen
        issue = 'noduplicate';
    end
    if sum(alist.I(e,p(pid),xcl) > 2.7 & alist.D(e,p(pid),xcl) > 1.5)
        xcell = C(e,p(pid),xcl);
    end
    if cellb < 0 || ~isnan(xcell)
        state = 'isdup';
    end
    
end

PC.SetFigure(DATA, DATA.tag.autoxyplot);
PC.PlotClusterXY(DATA, A{2});
PC.SetFigure(DATA,DATA.tag.autocelllist);
x = getappdata(gcf,'ShowDupBoxes');
for j = 1:length(x)
    if myhandle(x(j))
        delete(x(j));
    end
end
h(1) = PC.DrawBox(cells(1), cells(3), 'autolist','color','g');
h(2) = PC.DrawBox(cells(1), cells(2), 'autolist','color','g');
setappdata(gcf,'ShowDupBoxes',h);

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

yn = '';
if strcmp(issue,'duplicatea')
    ca = 1;
    cb = 2;
    icell = cella;
    probe = p(2);
else
    cb = 1;
    ca = 2;
    icell = cellb;
    probe = p(1);
end
if strcmp(state,'isdup')
    if truecell %both are duplicates of someone else
        s = sprintf('%s and %s both duplicates of %d(on p%.1f) Cl%salso duplicates.',str(A{ca}),str(A{cb}),-icell,truecell,sprintf('%d ',xcl));
    else
        s = sprintf('%s Already marked dup of Cell%d.  Cl%s also duplicates Matches two clusters in %s.',str(A{ca}),icell,sprintf('%d ',xcl),str(A{cb}));
    end
    yn = questdlg(s, 'Duplicate CheckPopup', 'OK','Next', 'OK');
elseif strcmp(state,'xcell')
    if truecell %both are duplicates of someone else
        s = sprintf('%s and %s both duplicates of %d(on p%.1f) Cl%s also duplicates.',str(A{ca}),str(A{cb}),-icell,truecell,sprintf('%d ',xcl));
    else
        s = sprintf('%s Already marked dup of Cell%d.  Cl%s (cell %d) also duplicates Matches two clusters in %s.',str(A{ca}),icell,sprintf('%d ',xcl),xcell,str(A{cb}));
    end
    yn = questdlg(s, 'Duplicate CheckPopup', 'Mark as Duplicate','OK','Next', 'OK');
    if strcmp(yn,'Mark as Duplicate')
        C (e,probe,xcl) = -(abs(icell));
        alist.CellLists{12} = C;
        setappdata(DATA.toplevel,'AutoCellList',alist);
        yn = 'Next';
    end
elseif strcmp(state,'isdup')
    s = sprintf('%s Not a cell.  Matches two clusters in %s.',str(A{1}),str(A{2}));
    yn = questdlg(s, 'Duplicate CheckPopup', 'OK','Next', 'OK');
end

if isempty(yn)
    response = yn;
else
    response = yn;
end

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

function MergePair(a,b)

line = get(a,'value');
strs = get(a,'string');
cells = sscanf(strs{line},'%d-%d');
DATA = GetDataFromFig(a);

function CheckGaps(a,b)

line = get(a,'value');
strs = get(a,'string');
pos = sscanf(strs{line},'Cell %d E%dP%d');
DATA = GetDataFromFig(a);


function FollowSquare(DATA, varargin)

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'square',5)
        j = j+1;
        square = varargin{j};
    end
    j = j+1;
end

FX = getappdata(DATA.toplevel,'AutoCellList');
Cs = FX.CellLists;
fprintf('  ');
for j = 1:size(square,1)
    fprintf('    E%dP%d.%d  ',square(j,:));
end
fprintf('\n');
cells = [];
for j = 1:length(Cs)
    if ~isempty(Cs{j})
        fprintf('%d:',j);
        c = celllist.get(Cs{j},square(1,:));
        fprintf('%s\n',sprintf('%9d   ',c));
        if c > 0
            cells(end+1) = c;
        end
    end
end

%X = FX.fixduplicates;
if isfield(FX,'fixboundary')
X = FX.fixboundary;
for k = 1:size(square,1)
for j = 1:length(X)
    id = find(X(j).elist == square(k,1) & X(j).probelist == square(k,2));
    if ~isempty(id)
        fprintf('Boundary fix %d\n',j);
    end
end
end
end

if ~isempty(cells)
    X = FX.cellq(cells(1));
    fprintf('Cell %d Started at E%dP%d.%d\n',cells(1),X.expt,X.probe,X.cl);
    for j = 1:length(X.matches)
        x = X.matches(j);
        if ~isempty(x.pos)
            fprintf('E%dP%d.%d %.3f I%.2f D%.2f\n',x.pos,x.xc,x.I,x.dropi)
        end
    end
end

function CheckAlignErrs(DATA, type, varargin)

FX = getappdata(DATA.toplevel,'AutoCellList')
CellList = PC.GetValue(DATA,'CellList');
x = FX.badtransitions;
for j = 1:length(x)
end
x = FX.alignerrs;
for j = 1:length(x)
    X = x(j);
    if sum(ismember(min(X.pos(1:2,1)):max(X.pos(1:2,1)),FX.shakelist))
        X.shake = 1;
    else
        X.shake = 0;
    end
    x(j).orderid = j;
    if x(j).aligned ==3 && length(X.newfit(2).matchdistance) > 1 && X.fixed <=0 && X.shake == 0
        matchd(j,:) = X.newfit(2).matchdistance(1:2);
        cld(j,:) = X.cldistance(1:2);
        cla(j,:) = X.cladistance(1:2);
        ncl(j) = length(X.cldistance);
        cells = CellList(X.pos(2,1),X.pos(2,2),X.duplicates(1:2));
        x(j).matchcells = squeeze(cells);
        cellid = find(cells>0);
        cells = cells(cellid);
        if length(cellid) > 1
            aligned(j) = 4;
        else
            aligned(j) = 3;
        end
    end
end
FX.alignerrs = x;
setappdata(DATA.toplevel,'AutoCellList',FX);
id = find(aligned ==4); %missing a split and its a cell
GetFigure('CompareDistances');
hold off;
plot(max(cld(id,:),[],2),max(matchd(id,:),[],2),'o');
set(gca,'xscale','log','yscale','log');
F = GetFigure('FixAutoList')
a = findobj(allchild(F),'flat','Tag','BoundaryList');
line = get(a,'value');
xid = find(max(matchd(id,:),[],2) < 3);
[distances,xid] = sort(max(matchd(id,:),[],2));
lx = get(a,'UserData');
allpos = cat(3,lx.pos);
good = getappdata(F,'goodlist');
oid = getappdata(F,'alignlist');
closeid = find(distances < 4);
if isempty(closeid)
    closeid = 1:length(distances)
end
for j = 1:length(closeid)
    X = x(id(xid(j)));
    marks{j} = SetMark(F,X,'BoundaryList');
end
%This finds fits where a split of C{1} gave good matches for two clusters
%ins C{2}, even if those two are not true duplicates of one in C{1}
setappdata(F, 'Markers',marks);
setappdata(F,'CurrentMarker',1)
X = x(id(xid(1)));
xpos = X.pos;
mid = find(allpos(1,1,:) == xpos(1,1) & allpos(2,1,:) == xpos(2,1) & allpos(1,2,:) == xpos(1,2));
mid = find([lx.orderid] == id(xid(j)));
if length(mid) ==1
    set(a,'value',mid);
    ShowBoundary(a, []);
    fprintf('? useful split');
end



