function X = PlotLFPpwr(X, varargin)
%
%Plot responses generated by PlotMLFP;
%
%[~,X] = PlotMLFP(LFP,'ftpwr','lines')
%   (If X is an LFP Expt, this line will be executed first, with default options)
% PlotLFPpwr(X,...)
% ...'print') prints a list of stimulus conditions to console
% ...'ids',[x]) only plots conditions in x
%               if x is a cell array, calculates the mean of ids in each
%               element of the array, then plots one line for each
% ...'showexamples' if X contains struct 'example', makes a plot for each of
%                  these. example(n).ids  sets what to plot
%...,'showexamples', example) 
%                            adds examples defined by the argument (which
%                            can combine across the cell arrays
% ...'showallallexmape' if X is a cell array, show examples for each element
%                       in additio to any defined by following arguemnt
% ...,'addexampe',X.plot); takes retured plot and adds it to the example
% list 

smw = 0;
colors = mycolors;
tight = 1;
printlabels = 0;
calcresp = 0;
coloroff = 0;
plotargs = {};
aplotargs = {}; %things to pass to set(gca)
pid = [];
explot = [];
holdon = 0;
showexamples = 0;
gtitle = '';
lfptag = 'LFPSpectrum';
calcfreq = [];
example = [];
examples = {};
showproperties = 0;
showallexamples = 0;
figurepos = [];
forcelabels = {};
linestyles = {};
legendpos = [];

if ischar(X)
    load(X);
    X = lp;
end
j =1;
while j <= length(varargin)
    if isfield(varargin{j},'ids')
        examples = {examples{:} varargin{j}};
    elseif strncmp(varargin{j},'addex',5)
        X = AddExample(X, varargin{j+1});
        return;
    elseif sum(strncmpi(varargin{j},{'calcresp' 'calcgamma'},5))
        calcresp = 1;
        if strncmpi(varargin{j},'calcgamma',5)
            calcfreq = [30 70];
        end
        X.plot.(varargin{j}) = varargin{j+1};
        if length(varargin) > j && ischar(varargin{j+1}) 
            j = j+1;
            if (isfield(X,varargin{j}) || iscell(X))
                explot = varargin{j};
                expts = 1;
            end
        elseif iscell(varargin{j+1});
            j = j+1;
            expts = varargin{j};
            explot = expts{1};
        end
    elseif strncmpi(varargin{j},'calcfreq',5)
        j = j+1;
        calcfreq = varargin{j};
    elseif strncmpi(varargin{j},'colors',3)
        j = j+1;
        colors = varargin{j};
        X.plot.colors = colors;
    elseif strncmpi(varargin{j},'coloroff',3)
        j = j+1;
        coloroff = varargin{j};
    elseif strncmpi(varargin{j},'hold',4)
        holdon = 1;
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        pid = varargin{j};
    elseif strncmpi(varargin{j},'labels',3)
        j = j+1;
        forcelabels = varargin{j};
    elseif strncmpi(varargin{j},'legendpos',7)
        j = j+1;
        legendpos = varargin{j};
    elseif strncmpi(varargin{j},'linestyles',7)
        j = j+1;
        linestyles = varargin{j};
    elseif strncmpi(varargin{j},'position',3)
        j = j+1;
        figurepos = varargin{j};
    elseif strncmpi(varargin{j},'print',3)
        printlabels = 1;
    elseif strncmpi(varargin{j},'title',5)
        j = j+1;
        gtitle = varargin{j};
    elseif sum(strcmpi(varargin{j},{'xlabel' 'ylabel' }))
        aplotargs = {aplotargs{:} varargin{j} varargin{j+1}};
        j = j+1;
    elseif strncmpi(varargin{j},'properties',6)
        showproperties = 1;
    elseif sum(strncmpi(varargin{j},{'showexamples' 'showallexamples'},6))
        showexamples = j; %record arg # so don't pass on this arg
        if strncmpi(varargin{j},'showall',6)
            showallexamples = 1;
        else
            showallexamples = 0;
        end
        if length(varargin) > j && isfield(varargin{j+1},'ids')
            j = j+1;
            example = varargin{j};
        end
    elseif strncmpi(varargin{j},'smooth',3)
        j = j+1;
        smw = varargin{j};
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        lfptag = varargin{j};
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end


if iscell(X)
    allf = {'n'};
    All = X{1};
    if isfield(All,'extypes')
        allf = All.extypes;
        for j = 2:length(X)
            allf = union(allf,[X{j}.extypes]);
        end
    else
        allf = {'n'};
    end
    allf = setdiff(allf,{''});
    idoffset = 0;
    All.example = [];
    for j = 1:length(X)
        for k = 1:length(allf)
            if isfield(X{j},allf{k})
                xf  = X{j}.(allf{k});
            elseif isfield(X{j},'Stimvals') && isfield(X{j}.Stimvals,allf{k})
                xf = ones(size(X{j}.n)) .* X{j}.Stimvals.(allf{k});
            else
                xf = zeros(size(X{j}.n));
            end
            if j == 1
                All.(allf{k}) =  xf;
            else
                All.(allf{k}) =  [All.(allf{k}) xf];
            end
        end
        if isfield(X{j},'example')
            n = length(All.example);
            for k = 1:length(X{j}.example)
                All.example = CatStruct(All.example,X{j}.example(k));
                All.example(k+n).ids = X{j}.example(k).ids + idoffset;
            end
        end
        idoffset = idoffset + size(X{j}.resps,3);
    end
    All.extypes = allf;
    for j = 2:length(X)
        All.labels = {All.labels{:} X{j}.labels{:}};
        All.n = [All.n  X{j}.n];
        All.resps = cat(3,All.resps, X{j}.resps);
    end
    if showexamples %first plot them all. Then add the new examples;
        res = PlotLFPpwr(All,varargin{1:showexamples-1});
    else
        res = PlotLFPpwr(All,varargin{:});
    end
    X = res;
    X.args = varargin;
    if showallexamples && isfield(X,'example')
        for j = 1:length(X.example);
            PlotExample(X,X.example(j),varargin{1:showexamples-1}); %args after 'showexample' only used for grouped data
        end
    end
    if isfield(example,'ids')
        X.example = CatStruct(X.example, example);
    if showexamples 
        for j = 1:length(example)
            PlotExample(X, example(j), varargin{1:showexamples-1});
        end        
    end
    end
    if isfield(res,'toplevel');
        AddGUI(res.toplevel,X);
    end
    return;
end


if ~isfield(X,'resps') && isfield(X,'Trials')  %not yet sent to PlotMLFP
    LFP = X;
%for manual Expts, just use the stimid    
    if strfind(LFP.Header.Options,'+exm')
        [~,X] = PlotMLFP(X,'ftpwr','lines','expts',{'stimid'});
    else
        [~,X] = PlotMLFP(X,'ftpwr','lines');
    end
end

if isempty(pid)
    pid = 1:size(X.resps,3);
elseif iscell(pid)
    Z = X;
    f = X.extypes;
    allid = [];
    
    for j = 1:length(pid)
        allid = [allid pid{j}(:)'];
        resps(:,1,j) = squeeze(mean(Z.resps(:,1,pid{j}),3));
        for k = 1:length(f)
            X.(f{k})(j) = mean(Z.(f{k})(pid{j}));
        end
    end
    props = ShowProps(X,allid);
    f = setdiff(fields(props),{'id'});
    for j = 1:length(pid)
        X.labels{j} = MakeLabel(X,f,pid{j});
    end
    X.resps = resps;
    pid = 1:size(X.resps,3);
end
if showproperties
    ShowProps(X,pid);
end

if isfield(X,'title') && isempty(gtitle)
    gtitle = X.title;
end
if ~isfield(X,'freqs')
    X.freqs = 1:size(X.resps,1);
end

[F, isnew] = GetFigure(lfptag,'front');
if isnew
    if isfield(X,'toplevel')  %a recursive call makeing new figure
        pos = get(X.toplevel,'position');
        pos(1) = pos(1)+pos(3);
        set(F,'position',pos);
    elseif ~isempty(figurepos)
        set(F,'position',figurepos);
    end
    X.toplevel = F;
end
X.pid = pid;
X.plot.tight = tight;
X.plot.ids = pid;

if calcresp
    for j = 1:size(X.resps,3)
        smoothed(:,j) = smooth(X.resps(:,1,j),5,'gauss');
    end
    if ~isempty(calcfreq)
        fid = find(X.freqs > calcfreq(1) & X.freqs < calcfreq(end));
        X.plot.calcfreq = calcfreq;
    else
        X.plot.calcfreq = calcfreq;
        rvar = std(smoothed,[],2)./mean(smoothed,2);
        [rmax, tmax] = max(rvar);
        x = prctile(rvar,[5 95]);
        crit = x(1) + diff(x) *0.15; %15% of max from baseline
        fid = find(rvar < crit);
        [rmax, tmax] = max(rvar);
        a = find(fid < tmax); %a(end) is last point below crit before peak
        b = find(fid > tmax);
        fid = fid(a(end))+1:fid(b(1));
    end
    baseline = min(smoothed(fid,:),[],2);
    for j = 1:size(X.resps,3)
        if isempty(calcfreq)
            X.gamma(j) = sum(X.resps(fid,:,j)-baseline)./sum(baseline);
        else
            X.gamma(j) = sum(X.resps(fid,:,j));
        end
    end
end
if smw
    for j = 1:size(X.resps,3)
        X.resps(:,1,j) = smooth(X.resps(:,1,j),smw,'gauss');
    end
    X.plot.smooth = smw;
end
if holdon == 0
    hold off;
end
baseline = min(smooth(squeeze(X.resps)',4,'gauss')',[],2);
for j  = 1:length(pid)
    if iscell(colors)
        c = colors{j};
    else
        c = colors(j);
    end
    h = plot(X.freqs,X.resps(:,1,pid(j)),'-','color',c,plotargs{:});
    if length(linestyles) >= j
        set(h,'linestyle',linestyles{j});
    end
    hold on;
end
legendlabels = X.labels(pid);
l = legend(legendlabels);
if ~isempty(legendpos)
    lpos = get(l,'position');
    lpos(1:length(legendpos)) = legendpos;
    set(l,'position',lpos);
end
%plot(baseline,'k');
set(gca,'yscale','log',aplotargs{:});
if tight
    axis('tight');
end

if ~isempty(gtitle)
    X.plot.title = gtitle;
    title(gtitle);
end

if printlabels
    for j = 1:length(X.labels)
        fprintf('%d: %s\n',j,X.labels{j});
    end
end


if ~isempty(explot)
    if strcmp(lfptag,'LFPSpectrum')
        [F, isnew] = GetFigure('GammaResp');
    else
        [F, isnew] = GetFigure([lfptag 'Gamma'],'front');
    end
    if isnew && isfield(X,'toplevel')
        pos = get(X.toplevel,'position');
        if pos(4)/2 > pos(2)
            pos(2) = pos(2) + pos(4);
        else
            pos(2) = pos(2) - pos(4);
        end
        set(F,'position',pos);
    end
    if holdon == 0
        hold off;
    else
        hold on;
    end
    if length(expts) > 1
        for j = 2:length(expts)
            ux = unique(X.(expts{j}));
            uvals{j} = ux;
            for k = 1:length(ux);
                ids{j-1}{k} = find(X.(expts{j}) == ux(k));
            end
            nv(j-1) = length(ux);
            cx{j-1} = 0;
        end
        nx = 0;
        allidx = [];
        for j = 1:prod(nv)
            [cx{:}] = ind2sub(nv,j);
            id = ids{1}{cx{1}};
            for k = 2:length(cx);
                id = intersect(id, ids{k}{cx{k}});
            end
            id = intersect(id,pid);
            id = setdiff(id, allidx);
            if ~isempty(id)
                nx = nx+1;
                allids{nx} = id;
                allidx = [allidx id];
                labels{nx} = '';
                for k = 1:length(cx)
                labels{nx} = sprintf('%s%s=%.2f ',labels{nx},expts{k+1},uvals{k+1}(cx{k}));
                end
            end
        end
    else
        allids{1} = pid;
        labels{1} = 'All';
    end
    for j = 1:length(forcelabels)
        if ~isempty(forcelabels{j})
            labels{j} = forcelabels{j};
        end
    end

    colors = mycolors;
    for j = 1:length(allids)
        c = j+coloroff;
        plot(X.(explot)(allids{j}),X.gamma(allids{j}),'o-','color',colors{c},...
            'markerfacecolor',colors{c});
        hold on;
    end
    legend(labels);
    xlabel(expts{1});
    title(sprintf('%.1f - %.1f Hz',X.freqs(fid(1)),X.freqs(fid(end))));
    GetFigure(lfptag,'noforce');
end
AddGUI(F, X);

if showexamples 
    if isempty(example) && isfield(X,'example') 
        example = X.example;
    end
    for j = 1:length(example)
        args = varargin(1:showexamples-1); %don't recurse
        PlotExample(X, example(j), args{:});
    end
    return;
end

function X = AddExample(X, plot, varargin)
% create an example from current plot
    cfields = {'ids' 'calcresp' 'calcfreq' 'calcgamma' 'colors' 'labels' 'linestyles' 'title'};
if isfield(X,'example')
    X.example(end+1).title = 'New';
else
    X.example(1).title = 'New';
end

for k = 1:length(cfields)
    if isfield(plot,cfields{k}) && ~isempty(plot.(cfields{k}))
        X.example(end).(cfields{k}) = plot.(cfields{k});
    end
end

function PlotExample(X, example, varargin)

args = varargin;
    cfields = {'calcresp' 'calcfreq' 'calcgamma' 'colors' 'labels' 'linestyles'};
    plotfields = {'smooth' 'calcfreq'};
if isfield(example,'title') && ~isempty(example.title)
    tag = example.title;
    tstr = example.title;
else
    fs = findobj('type','figure');
    tags = get(fs,'tag');
    n = sum(strncmp('Spectrum',tags,7));
    tag = sprintf('Spectrum%d',n);
    tstr = [];
end
for k = 1:length(plotfields)
    if isfield(X.plot,plotfields{k})
        args = {args{:} plotfields{k} X.plot.(plotfields{k})};
    end
end
for k = 1:length(cfields)
    if isfield(example,cfields{k}) && ~isempty(example.(cfields{k}))
        args = {args{:} cfields{k} example.(cfields{k})};
    end
end
PlotLFPpwr(X, args{:}, 'ids',example.ids,'tag',tag);
if ~isempty(tstr)
    title(tstr);
end
if isfield(example,'Comment')
    for j = 1:length(example.Comment)
        fprintf('%s\n',example.Comment{j});
    end
end


function AddGUI(F, X)
hm = findobj(F,'tag','PlotMenu');
if isempty(hm)
    hm  = uimenu(F,'label','Plot','Tag','PlotMenu');
else
    delete(get(hm,'children'));
end
X.selected = zeros(1,length(X.labels));
if ~isfield(X,'pid')
    X.pid = 1:length(X.selected);
end
X.selected(X.pid) = 1;
setappdata(F,'LFPresps',X);

uimenu(hm,'Label','print','callback',{@PlotSpectrum, 'print', 1});
for j = 1:length(X.labels)
    h = uimenu(hm,'Label',X.labels{j},'callback',{@PlotSpectrum, 'select', j});
    if X.selected(j)
        set(h,'checked','on');
    end
end

hm = findobj(F,'tag','ExampleMenu');
if isempty(hm)
    hm  = uimenu(F,'label','Examples','Tag','ExampleMenu');
else
    delete(get(hm,'children'));
end
if isfield(X,'example')
for j = 1:length(X.example)
    str = X.example(j).title;
    if isempty(str)
        str = sprintf('Set%d',j);
    end
    h = uimenu(hm,'Label',str,'callback',{@PlotSpectrum, 'example', j});
end
end


function PlotSpectrum(a,b,type,ind)

F = GetFigure(a);
X = getappdata(F,'LFPresps');

onoff = {'off' 'on'};
if strcmp(type,'select')
X.selected(ind) = ~X.selected(ind);
set(a,'checked',onoff{1+X.selected(ind)});
setappdata(F,'LFPresps',X);
hold off;
pid = find(X.selected);
baseline = min(smooth(squeeze(X.resps)',4,'gauss')',[],2);
colors = mycolors;
for j  = 1:length(pid)
    c = colors{j};
    plot(X.freqs,X.resps(:,1,pid(j)),'-','color',c);
    hold on;
end
legend(X.labels(pid));
set(gca,'yscale','log');
if X.plot.tight
    axis('tight');
end
elseif strcmp(type,'example')
    calcfreq = [];
    args = {};
    E = X.example(ind);
    PlotExample(X, E);
elseif strcmp(type,'print')
    for j = 1:length(X.labels)
        fprintf('%d: %s\n',j,X.labels{j});
    end
end

function props = ShowProps(X, id)

props.id = id;
for j = 1:length(X.extypes)
    f = X.extypes{j};
    if length(unique(X.(f)(id))) > 1
        props.(f) = X.(f)(id);
        fprintf('%s %s\n',f,sprintf(' %.2f',X.(f)(id)));
    end
end

function lbl = MakeLabel(X, f, id)

lbl = '';

if iscell(X)
    for k = 1:length(f)
        lbl = sprintf('%s%s=%.2f ',labels{nx},expts{k},uvals{k}(f{k}));
    end
else
    for j = 1:length(f)
        lbl = sprintf('%s%s=%.2f ',lbl,f{j},mean(X.(f{j})(id)));
    end
end
