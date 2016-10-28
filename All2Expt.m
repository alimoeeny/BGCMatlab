function Expt = All2Expt(AllE, cell, varargin)
%Expt = All2Expt(AllExpt, cell, varargin)
%makes a single Expt struct from a multicell AllExpt struct
%returns Expt for cellnumber == cell (not necessarily elemenn# == cell)
%All2Expt(AllE, n,'all') returns data for nth element, regardless of
%
%All2Expt(AllE, n,'withmu') For cases where a cell is only present for a
%subset of experiments, this combines the MU expts for the same probe with
%the SU. Expt.Header.iscell shows which blocks are SU.
%(In these cases, combine sets the rate in MU blocks to match the mean in
%the SU blocks). 
%
%Expts = All2Expt(AllExpt) returns a cell array of Expts
%
mu = 0;
addmu = 0;
makeallcells = 0;
j = 1;

narg = nargin;
if nargin > 1 && ischar(cell)
    narg = 1;
    varargin = {cell varargin{:}};
end

while j <= length(varargin)
    if strncmpi(varargin{j},'allcells',8)
        mu= 1;
    elseif strncmpi(varargin{j},'all',2)
        mu= 2;
    elseif strncmpi(varargin{j},'mu',2)
        mu= 1;
    elseif strncmpi(varargin{j},'probe',2)
        findprobe =1;
    elseif strncmpi(varargin{j},'withmu',2)
        addmu = 1;
    end
    j = j+1;
end

Expt = AllE.Expt;

cells = [AllE.Header.cellnumber];
cells = cells((cells > 0));
probes = [AllE.Header.probe];
if addmu && isfield(AllE.Header,'muforcell')
    for j = 1:length(AllE.Header)
        [a,b] = Counts(AllE.Header(j).muforcell,'descend');
        if ~isempty(a)
            mucell(j) = b(1);
        else
            mucell(j) = 0;
        end
    end
end

if narg < 2 || ischar(cell)
    if mu ==2
        id = 1:length(AllE.Header);
    else
        id = 1:length(cells);
        makeallcells = 1;
    end
elseif mu == 0
    id = find(ismember(cells, cell));
    if length(cell) == 1 && length(id) > 1
        cprintf('red','Two Expts matching Cell %d in %s\n',cell,GetName(AllE,'loadname'));
        id = id(1);
    end
elseif mu == 2
    id = cell;
elseif findprobe
    id = find(probes == cell);
    id = id(1); %is cell and MU, use cell
else
    id = find(probes == cell & cells == 0);
end

if isempty(id)
    [Expt.Trials.Spikes] = deal([]);
    if makeallcells
        Expts{1} = Expt;
        Expt = Expts;
    end
    return;
end
if length(id) > 1 || makeallcells %make multiple expts as cell array
    for j = 1:length(id)
        Expts{j} = All2Expt(AllE,id(j), 'all');
    end
    Expt = Expts;
    return;
elseif mu == 2 %asked for all cells but only one here. Still return a cell array
%    Expts{1} = Expt;
%    Expt = Expts;
end
    
%fields that are not stored in each Header, but depend on the blocks
%included
blockfields = {'BlockStart','BlockStartid' 'depths' 'BlockTimeOffset'};
if isfield(AllE.Header,'Combineids')
    bid = find(ismember(Expt.Header.Combineids,AllE.Header(id).Combineids));
    for j = 1:length(blockfields)
        if isfield(AllE.Expt.Header,blockfields{j})
            Expt.Header.(blockfields{j}) = AllE.Expt.Header.(blockfields{j})(bid);
        end
    end
end
f = setdiff(fields(AllE.Header),blockfields);
for j = 1:length(f)
    Expt.Header.(f{j}) = AllE.Header(id).(f{j});
end

if ~isempty(Expt.Header.Clusters) && iscell(Expt.Header.Clusters{1})
    for j = 1:length(Expt.Header.Clusters)
        p = round(Expt.Header.probe); %temp kludge  Need to figure out probe by block
        Expt.Header.Clusters{j} = Expt.Header.Clusters{j}{p};
    end
    Expt.Header.DataType =  'Spike2Swatches';
end

cellno = AllE.Header(id).cellnumber;
p = AllE.Header(id).probe;
if isfield(AllE.Spikes{id},'Trial')
    [~, uset, usespk] = intersect([Expt.Trials.Trial],AllE.Spikes{id}.Trial);
else
    [~, uset, usespk] = intersect([Expt.Trials.id],AllE.Spikes{id}.trialid);
end
if length(uset) ~= length(AllE.Spikes{id}.trialid)
    missing = setdiff(AllE.Spikes{id}.trialid,[Expt.Trials.id]);
    emissing = setdiff([Expt.Trials.id], AllE.Spikes{id}.trialid);

   Expt = AddError(Expt,'%s Spikes %d Cell%d P%d Trial Id length mismatch Missing %s',GetName(Expt),id,cellno,p,sprintf(' %d',missing));
   AllE.Spikes{id}.Spikes = AllE.Spikes{id}.Spikes(usespk); 
end

Expt.Trials = Expt.Trials(uset);
Expt.trialsused = uset;
counts = zeros(1,length(AllE.Spikes{id}.Spikes));
for j = 1:length(AllE.Spikes{id}.Spikes)
    Expt.Trials(j).Spikes = double(AllE.Spikes{id}.Spikes{j});
    Expt.Trials(j).OSpikes = double(AllE.Spikes{id}.OSpikes{j});
    Expt.Trials(j).Ocodes = double(AllE.Spikes{id}.Ocodes{j});
    count = sum(Expt.Trials(j).Spikes > 500 & Expt.Trials(j).Spikes < Expt.Trials(j).dur+500);
    Expt.Trials(j).count = count;
    counts (j) = count;
end

if isfield(AllE.Header,'excludeids') && ~isempty(AllE.Header(id).excludeids)
    %NB don't use setdiff in case there are repeats in Trials.id
    gid = find(~ismember([Expt.Trials.id],AllE.Header(id).excludeids));
    Expt.Trials = Expt.Trials(gid);
end
if ~isfield(Expt.Header,'nspk')
    Expt.Header.nspk = sum(counts);
end

if addmu
    id = find(mucell == cell)
    if ~isempty(id)
        muExpt = All2Expt(AllE, id(1),'all');
        Expt = CombineMU(Expt,muExpt);
    end
end



function Expt = CombineMU(Ea, Eb)

Expt = Ea;
T = cat(2,Ea.Trials,Eb.Trials);
[a,b] = sort([T.Trial]);
T = T(b);
Expt.Trials = T;
H = Ea.Header;
f = {'BlockStart' 'suffixes' 'Combineids' 'Clusters' 'excludeids' 'dips' 'dropi' 'BlockStartid' 'depths'};

for j = 1:length(f)
    if isfield(H,f{j}) && isfield(Eb.Header,f{j})
        if strcmp(f{j}, 'dips')
            bs = cat(1,Ea.Header.(f{j}),Eb.Header.(f{j}));
            bs = bs(sortid,:);
        else
            bs = cat(2,Ea.Header.(f{j}),Eb.Header.(f{j}));
        end
        if j == 1  %get order
            [a, sortid] = sort(bs);
        end
        if sum(strcmp(f{j},{'excludeids' 'dips'})) % can't use sort list as # may vary
            H.(f{j}) = bs;
        else
            H.(f{j}) = bs(sortid);
        end
    end
end
H.Start = min([Ea.Header.Start Eb.Header.Start]);
H.End = max([Ea.Header.End Eb.Header.End]);
H.iscell = ismember(H.BlockStart,Ea.Header.BlockStart);
Expt.Header = H;
