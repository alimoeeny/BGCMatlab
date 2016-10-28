function [result] = PlotRateSequence(Expts, varargin)
%PlotRateSequence(Expts, ...) plots the rates for each trial in a set of expts
%(different blocks for one cell)
%PlotRateSequence(AllExpt) Plots sequences for each element in AllExpt
%PlotRateSequence(AllExpt,'cells') Plots sequences for each cell in AllExpt
% see also CheckExptRates

tn = [];
counts = [];
AllTrials = [];
AllIds = [];
AllTid = [];
result = [];
AllBlocks = [];
xstyle = 'trial';
normalize = 0;
   offset = 0;
   color = 'b';
   callback = {};
plotcells = 0;
checkcounts = 0;
plottype = 'rate';
state.tag = 'MeanRates';
smoothw = 0;

   j = 1;
   while j <= length(varargin) 
       if strncmpi(varargin{j},'bytime',5)
           xstyle = 'time';
       elseif strncmpi(varargin{j},'cells',5) 
%Expt expt{} is a different cell
% or if given and AllExpt struct, only plot cells
           plotcells = 1;
       elseif strncmpi(varargin{j},'blockcount',8)
           plottype = 'blockrate';
       elseif strncmpi(varargin{j},'check',5)
           checkcounts = 1;
       elseif strncmpi(varargin{j},'callback',8)
           j = j+1; 
           callback = varargin{j};
       elseif strncmpi(varargin{j},'color',5)
           j = j+1; 
           color = varargin{j};
       elseif strncmpi(varargin{j},'normalize',5)
           normalize = 1;
       elseif strncmpi(varargin{j},'offset',5)
           j = j+1;
           offset = varargin{j};
       elseif strncmpi(varargin{j},'smooth',5)
           j = j+1;
           smoothw = varargin{j};
       end
       j = j+1;
   end
   
%?? this causes infinite recusion with iscell() below
%why was it here
if isstruct(Expts) && isfield(Expts,'Trials') && 0
    Expt = Expts;
    clear Expts;
    Expts{1} = Expt;
end
if isstruct(Expts) && isfield(Expts,'allrates')
    GetFigure(state.tag);
    result = CheckCounts(Expts, state, varargin{:});
    return;
end
if iscellstr(Expts)
    names = Expts;
    for j = 1:length(names)
        Expt = LoadExpt(names{j});
        if isempty(Expt)
            result{j}.name = names{j};
            AddError(result{j},'Empty Expt in %s',names{j});
        else
            result{j} = PlotRateSequence(Expt,varargin{:});
        end
    end
    return;
elseif iscell(Expts) && ~isempty(CellToMat(Expts,'maxtrial')) %array of results
    for j = 1:length(Expts)
        result{j} = PlotRateSequence(Expts{j},varargin{:});
    end
    return;
elseif ischar(Expts)
    
    name = Expts;
    Expt = LoadExpt(name);
    result = PlotRateSequence(Expt,varargin{:});
    result.Expt = Expt;
return;
end


if isstruct(Expts) && isfield(Expts,'Spikes');
    A = Expts;
    clear Expts;
    if plotcells
        Expts = All2Expt(A);
    else
        Expts = All2Expt(A,'all');
    end
    filename = A.Expt.Header.loadname;
    plotcells = 1;
    result.maxtrial = max([A.Expt.Trials.Trial]);
end

cellnumber = [];
if plotcells
    colors = mycolors;
    y = 0;
    hold off;
    for j = 1:length(Expts)
        if isfield(Expts{j},'plotres')
            E = Expts{j}.plotres.Data;
        else
            E = Expts{j};
        end
        if ~isempty(E.Trials)
        if normalize
            x{j} = PlotRateSequence(E,'color',colors{j}, 'offset', (j-1),'normalize');
            line(get(gca,'xlim'),[j-1 j-1],'color',colors{j});
        else
            x{j} = PlotRateSequence(E,'color',colors{j}, 'offset', y);
            line(get(gca,'xlim'),[y y],'color',colors{j});
            y = y+max([E.Trials.count]);
        end
        if isfield(E,'trialsused')
            x{j}.uset = E.trialsused;
        else
            x{j}.uset = 1:length(x{j}.rates);
        end
        x{j}.ids = [E.Trials.id];
        allids{j} = x{j}.times(:);
        cellnumber(j) = E.Header.cellnumber;
        hold on;
        end
    end
    if checkcounts
        ntrials = length(unique(cat(1,allids{:})));
        alltrials = CellToMat(x,'times'); %trial #, usually
        triallist = unique(alltrials);
        triallist = triallist(triallist > 0);
        allrates = ones(length(x),ntrials).*NaN;
        for j = 1:length(x)
            allrates(j,x{j}.uset) = x{j}.rates - x{j}.offset;
            ids(x{j}.uset) = x{j}.ids;
        end
        result.allrates = allrates;

        a = sum(allrates ==0);
        id = find(a > length(x)/4);
        if ~isempty(id)
            result = AddError(result,'%s: Zero rates on %d-%d probes Trials %s\n',GetEval(E,'shortname'),min(a(id)),max(a),sprintf(' %d',ids(id)));
        end
        result.trials = triallist;
        result.ids = ids;
        result.missingid = ids(id);
        result.missing = id;
        result.name = GetName(E,'path');
        result.checktime = now;
        triallist = triallist(triallist > 0);
        result.cellnumber = cellnumber;
        for j = 1:length(triallist)
            id = find(alltrials == triallist(j));
            meanrate(j) = mean(allrates(id));
            meanvar(j) = var(allrates(id));
        end
        GetFigure(state.tag);
        hold off;
        plot(triallist,meanrate,'k-');
        hold on;
        plot(triallist,meanvar,'r-');
    end
    AddPlotMenu(gcf);
    return;
end

celldef = [];
%This plots a list of Expts for one cell.
   for j = 1:length(Expts)
       if iscell(Expts)
           Expt = Expts{j};
       else
           Expt = Expts(j);
       end
       E = Expt;
       if isfield(E.Header,'cellnumber')
            cellnumber(j) = E.Header.cellnumber;
       else
           cellnumber(j) = 0;
       end

    dur = mean([Expt.Trials.dur])./10000;
    rates = [Expt.Trials.count]./dur;
    counts = [counts rates];
    celldef(end+1:length(counts)) = cellnumber(j);
    result.meanrates(j) = mean(rates);
    tn = [tn Expt.Trials.Trial];
    ends(j) = tn(end);
        AllIds = cat(2,AllIds,[E.Trials.id]);
        if strcmp(xstyle,'time')
            if isfield(E.Header,'timeoffset')
                toff = E.Header.timeoffset;
                if isfield(E.Header,'timeadjust')
                    toff = E.Header.timeoffset-E.Header.timeadjust;
                end
            else
                toff = 0;
            end
            AllTrials = cat(2,AllTrials,toff+([E.Trials.TrialStart]./10000));
            if isempty(E.Trials)
                AllBlocks(j) =  E.Header.trange(1)+toff;
            else
                AllBlocks(j) =  E.Trials(1).TrialStart./10000 +toff;
            end
            AllTid = cat(2,AllTid,[E.Trials.Trial]);
        else
            AllTrials = cat(2,AllTrials,[E.Trials.Trial]);
            AllTid = AllTrials;
            AllBlocks(j) =  E.Trials(1).Trial;
            if isfield(E.Header,'cellnumber')
                cellnumber(j) = E.Header.cellnumber;
            end
        end
   end
   if length(Expts) == 1 && isfield(Expt.Header,'BlockStart')
       AllBlocks = Expt.Header.BlockStart;
   end
   result.times = AllTrials;
%if there are duplicate Trials ?? why do we delete some?   
%should not happen in files where Trials.Trial has been fixed
    [AllT, id] = unique(AllTrials);
    [uniqueid, uid] = unique(AllIds);
    if length(id) < length(result.times)
        if length(uid) <= length(id)
            result = AddError(result,'-show','PlotRateSequence: %d Duplicate Trials',length(result.times)-length(id));
        elseif length(uid) < length(result.times)
            result = AddError(result,'-show','PlotRateSequence: %d Duplicate Ids',length(result.times)-length(uid));
            id = uid;
        else
            result = AddError(result,'-show','PlotRateSequence: %d Duplicate Trials but Ids unique',length(result.times)-length(id));            
            id = 1:length(AllIds);
        end
    end
    AllTrials = AllTrials(id);
    gid = find(AllTrials > 0);
    AllTrials = AllTrials(gid);
    result.times = AllTrials;

    AllIds = AllIds(id(gid));
    AllTid = AllTid(id(gid));
    counts = counts(id(gid));

    
    dy = 0;
    if normalize == 1
        scale = 1./mean(counts);
        dy = offset;
    elseif normalize == 2
        dy= mean(counts) .* offset;
        scale = 1;
    else
        dy = offset;
        scale = 1;
    end
    if ~isempty(AllBlocks)
        blockid = zeros(1,length(Expt.Trials));
        for j = 1:length(AllBlocks)-1
            t = AllBlocks(j):AllBlocks(j+1);
            [~, tid] = intersect([Expt.Trials.Trial],t);
            blockid(tid) = j;
        end 
        t = AllBlocks(end):Expt.Trials(end).Trial;
        [~, tid] = intersect([Expt.Trials.Trial],t);
        blockid(tid) = length(AllBlocks);
        blk = unique(blockid);
        for j = 1:length(blk)
            tid = find(blockid == blk(j));
            result.blockrates(j) = mean(counts(tid));
        end
    end
    result.meanrates = (result.meanrates .* scale) + dy;
    result.rates = (scale.*counts)+dy;
    result.offset = dy;
    result.scale = scale;
    result.cellnumbers = cellnumber;
    result.meancount = mean(counts);
    result.ids = AllIds;
    result.Trials = AllTid;
    result.options.action = 'print';
    
    if isempty(callback)
        if strcmp(plottype,'blockrate')
            h = plot(AllBlocks,result.blockrates,'o-','color',color);
        else
            h  = plot(AllTrials, (scale.*counts)+dy,'o','color',color,'buttondownfcn',@HitRatePlot);
        end
        setappdata(gcf,'RateSeq',result);
        setappdata(gcf,'RateExpts',Expts);
    else
        h  = plot(AllTrials, (scale.*counts)+dy,'o','color',color,'buttondownfcn',callback);
    end
    id = find(celldef ==0);
    if ~isempty(id) && sum(celldef) > 0 %includes MU blocks
        hold on;
        h  = plot(AllTrials(id), (scale.*counts(id))+dy,'o','color',color,'markerfacecolor',color);
        if ~isempty(callback)
            set(h,'buttondownfcn',callback);
        else
            set(h,'buttondownfcn',@HitRatePlot);
        end
    end
    if smoothw
        hold on;
        h  = plot(AllTrials, smooth((scale.*counts)+dy,smoothw),'-','color',color);
    end
    if ~strcmp(plottype,'blockrate')
    for j = 1:length(AllBlocks)        
        hold on;
        plot([AllBlocks(j) AllBlocks(j)], get(gca,'ylim'),'--');
    end
    end
    if length(Expts) == 1
        title(sprintf('E%d %s',GetExptNumber(Expts),GetName(Expts,'withcell','withexpt')));
    end
    result.AllBlocks = AllBlocks;
    result.handle = h;
    AddPlotMenu(gcf);

function x = CheckCounts(x, state, varargin)
    
    cellsonly = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cellsonly',4)
            cellsonly = 1;
        end
        j = j+1;
    end
    showsd = 0;
    a = sum(x.allrates ==0);
    id = find(a > size(x.allrates,1)/4);
    if ~isempty(id)
        fprintf('%s: Zero rates on %d-%d probes Columns %s\n',GetName(x),min(a(id)),max(a),sprintf(' %d',x.ids(id)));
    end
    %x.missingid = ids(id);
    x.missing = id;
    x.checktime = now;
    if cellsonly && isfield(x,'cellnumber')
        cid = find(x.cellnumber > 0);
    else
        cid = 1:size(x.allrates,1);
    end
    triallist = x.trials;
    for j = 1:size(x.allrates,2)
        meanrate(j) = nanmean(x.allrates(:,j));
        meanvar(j) = nanvar(x.allrates(:,j));
    end
    GetFigure(state.tag);
    hold off;
    id = find(~isnan(meanrate));
    plot(triallist(id),meanrate(id),'k-');
    if showsd
        hold on;
        plot(triallist(id),meanvar(id),'r-');
    end
    x.zrange = minmax(zscore(meanrate(id)));

    
    function HitPoint(a,b,t)
        
        fprintf('Trial %d\n',t);
        
function [pos, id] = FindNearestPoint(line, x, y)
    ax = gca;
    pos = get(ax,'currentpoint');
    xl = get(ax,'xlim');
    yl = get(ax,'ylim');
    r = (pos(1,1) -x)./diff(xl) + i .* (pos(1,2) -y)./diff(yl);
    [a,b] = min(abs(r));
    pos = [x(b) y(b)];
    id = b;

    
function KeepSpks(F, Spks)    
%KeepSpks(F, Spks) store Spikes in appdata
AllSpks = getappdata(F,'AllSpks');
for k = 1:length(Spks)
    new = 1;
    for j = 1:length(AllSpks)
        if Spks{k}.exptno == AllSpks{j}.exptno
            new = 0;
        end
    end
    if new
        AllSpks{end+1} = Spks{k};
        setappdata(F,'AllSpks',AllSpks);
    end
end


function HitRatePlot(data,b)
    R = getappdata(gcf,'RateSeq');
    E = getappdata(gcf,'RateExpts');
    F = gcf;
    AllSpks = getappdata(F,'AllSpks');
    [x, id] = FindNearestPoint(data, R.times, R.rates);
    fprintf('Trial %d id%d (#%d)\n',R.Trials(id),R.ids(id),id);
    if strcmp(R.options.action,'quickspks')
        Spks = expt.PlotSpikes(E,'quickspks','Trial',R.Trials(id),AllSpks);
    elseif strcmp(R.options.action,'trialspks')
        Spks = expt.PlotSpikes(E,'selecttrial',R.Trials(id),AllSpks);
    elseif strcmp(R.options.action,'spool')
        tid = expt.FindBlocks(E,R.Trials(id),'trials'); %list of trials in block
        Spks = expt.PlotSpikes(E,'selecttrials',tid,AllSpks);
    else
        Spks = [];
    end
    KeepSpks(F, Spks);
        
    
function MenuHit(src, b, type, fcn)    
    
    F = GetFigure(src);
    R = getappdata(gcf,'RateSeq');
    if strcmp(type,'action')
        SetMenuCheck(src,'exclusive');
        R.options.action = fcn;
        setappdata(F,'RateSeq',R);
    end
    
function AddPlotMenu(F, varargin)
    x = findobj(allchild(F),'flat','type','uimenu','tag','RatePlotMenu');
    
    if ~isempty(x)
        return;
    end
    hm = uimenu(F,'label','Options','tag','RatePlotMenu');
    sm = uimenu(hm,'label','Callback Action');
    uimenu(sm,'label','Just Print Trial #','callback',{@MenuHit, 'action', 'trial'});
    uimenu(sm,'label','Spool Block','callback',{@MenuHit, 'action', 'spool'});
    uimenu(sm,'label','Quick Spikes for Block','callback',{@MenuHit, 'action', 'quickspks'});
    uimenu(sm,'label','Trial Spikes','callback',{@MenuHit, 'action', 'trialspks'});
    
    
    
    
    