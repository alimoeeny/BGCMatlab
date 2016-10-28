function [names, details] = ListExpts(E, varargin)
%[names, details] = ListExpts(E, ...)   lists names of Expts in  E
% E can be a cell array of expts, or a directory name
% If its  directory name, only that directorly is listed
% ListExpts(E, 'depth',1) lists all directoris one layer down 
% For Directories, ListExpts Lookd in .idx files for experiments performed.
% ListExpts can search its own output to match specific expt names, and 
% then to search for combined Expt files with a given suffix
%
%ListExpts(dir,'check')  checks to see if any files have been updated
%since last listing was made
%
%ListExpts(...,'show','xx')
%      shows the value for field 'xx' in each expt. If this varies, the
%      mean value is printed
%
% [names, details] = ListExpts(dir);
% list = ListExpts(name, 'rds.dxXce', details);
% list = ListExpts(details, 'rds.dxXce');
%      equivalent. return list of expts that match
%         returns a structure with a list of directoies that contain this
%         expt
% exlist = ListExpts(list,'AC');
%    Counts the numbe rof matching files in each folder, and the number of
%    cell files
% see also FindExpts - easier for finding expts of one type

depth = 0;
psych = 0;
showvals = {};
j = 1;
argson = {};
while j <= length(varargin)
    if strncmp(varargin{j},'depth',5)
        j = j+1;
        depth = varargin{j};
    elseif strncmp(varargin{j},'psych',5)
        psych = 1;
    elseif strncmp(varargin{j},'show',5)
        j = j+1;
        showvals{end+1} = varargin{j};
    end
    j = j + 1;
end

if iscellstr(E) %list of names
    if isdir(E{1})
        for j = 1:length(E)
            [names{j}, details{j}] = ListExpts(E{j}, varargin{:});
        end
        return;
    else
        [names, details] = FindExpts(E, varargin{:});
    end
elseif psych
    [names,details] = ListPsychDir(E, 0);
elseif iscell(E)
    if iscellstr(E{1}) %E is set of cellstrs, with names for each dir
        [names, details] = FindExpts(E, varargin{:});
    elseif isfield(E{1},'dirpath')
        [names, details] = FindExpts(E, varargin{:});
    else
        [names, details] = ListExptCells(E, varargin{:});
    end
elseif isstruct(E)  % a details struct with a list of dirs containing expts. 
    if isfield(E,'dirpath') && iscell(E.dirpath) 
        names = CheckForExpts(E, varargin{:}); %first arg should be suffix to check
    elseif isfield(E,'dirpath') %just one
        if ~isfield(E,'matchnames')
            dirpath = E.dirpath;
            names = E.names;
            E = rmfields(E,{'dirpath' 'names'});    
            E.names{1} = names;
            E.dirpath{1} = dirpath;
            E.matchnames{1} = unique(names);
            for k = 1:length(E.matchnames{1})
                id = find(strcmp(E.matchnames{1}{k},names));
                E.trialcounts(1,k) = sum(E.ntrials(id));
            end
        end
        names = CheckForExpts(E, varargin{:}); %first arg should be suffix to check
    end
elseif isdir(E)
    if depth > 0
    [names, details] = ListExptDirs(E, depth, varargin{:});
    else
    [names, details] = ListExptDir(E, varargin{:});
    end
elseif ischar(E)
    [names, details] = ListExptDirs(E, depth, varargin{:});
end


function E = CheckForExptsByType(E, type, varargin)

if ~isfield(E,'params')
    return;
end

for j = 1:length(E.params)
    if sum(strcmp(type,E.params{j}))
        suffixes{j} = MakeName(E,j);
    else
        suffixes{j} = '';
    end
end
suffixes = setdiff(unique(suffixes),'');
for j = 1:length(suffixes)
    Es{j} = CheckForExpts(E, suffixes{j});
end

function DiskFiles = CheckForExpts(E, suffix, varargin)
%Finds disk files with combined data
extype = '';
condense = 1; %default is not to keep raw data in fit results
if strcmp(suffix,'type')
    extype = varargin{1};
end
strargs = cell2cellstr(varargin);
if sum(strcmp('fit',strargs))
    fitcells = 1;
else
    fitcells = 0;
end

for e = 1:length(E)
    if iscell(E(e).dirpath)
        dirpath = E(e).dirpath;
    else
        dirpath{1} = E(e).dirpath;
    end
    if ~isempty(extype)
        E(e) =CheckForExptsByType(E(e),extype);
        dirpath = '';
    end
    for j = 1:length(dirpath)
        if sum(CellToMat(regexp(E(e).matchnames{j},suffix)));
            good(j) = 1;
        end
    end
    good = find(good);
    D.pattern = suffix;
    for j = 1:length(good);
        did = good(j);
        if sum(strcmp('exact',strargs))
            s  = [dirpath{good(j)} '/*' suffix '.*'];
        else
            s  = [dirpath{good(j)} '/*' suffix '*'];
        end
        d = dir(s);
        xsuff = expt.Name2Suffix(suffix);
        if ~isempty(xsuff)
            if sum(strcmp('exact',strargs))
                s  = [dirpath{good(j)} '/*' xsuff '.*'];
            else
                s  = [dirpath{good(j)} '/*' xsuff '*'];
            end
            db = dir(s);
            d = cat(1,d,db);
        end
        D.files{j} = {d.name};
        D.nfiles(j) = length(d);
        D.dirpath{j} = dirpath{did};
        blkid = find(CellToMat(regexp(E(e).names{did},suffix)));
        if sum(strcmp('exact',strargs))
            blkid = find(CellToMat(regexp(E(e).names{did},[suffix '$'])));
            suffid = find(CellToMat(regexp(E(e).matchnames{did},[suffix '$'])));
            D.nblock(j)= length(suffid);
        else
            suffid = find(CellToMat(regexp(E(e).matchnames{did},suffix)));
        end
        D.nblock(j)= length(blkid);
        cells = [];
        xcells = [];
        cellsfile = '';
        for k = 1:length(d)
            id  =findstr(d(k).name,'Cells');
            if ~isempty(id)
                fprintf('Checking %s\n',d(k).name);
                cellsfile = d(k).name;
                AllExpt = LoadExpt([dirpath{did} '/' d(k).name]);
                for c = 1:length(AllExpt.Header)
                    H = AllExpt.Header(c);
                    if H.cellnumber > 0
                        cells = [cells H.cellnumber];
                        D.celltrials{j}(H.cellnumber) = length(AllExpt.Spikes{c}.trialid);
                        D.nspikes{j}(H.cellnumber) = length(cat(1,AllExpt.Spikes{c}.Spikes{:}));
                        D.durations{j}(H.cellnumber) = expt.value(AllExpt.Expt,'duration');
                        if fitcells
                            res = PlotExpt(All2Expt(AllExpt,H.cellnumber),'frbox','noplot','rcnmin',10,'splitsdf');
                            efit = FitExpt(res);
                            if isfield(efit,'RespVar');
                                D.respvar{j}(H.cellnumber) = efit.RespVar;
                            elseif isfield(efit,'di')
                                D.di{j}(H.cellnumber) = efit.di;
                            end
                            D.fits{j}{H.cellnumber} = efit;
                            if condense >= 1
                                D.plots{j}{H.cellnumber} = rmfields(res,{'Data'});
                            else
                                D.plots{j}{H.cellnumber} = res;
                            end
                        end
                        D.cellfiles{j}{H.cellnumber} = GetName(AllExpt,'loadname');
                    end
                end
            end            
        end
        ncells = length(cells);
        filename = cellsfile;
        for k = 1:length(d)
            cellnumber = 0;
            id  =findstr(d(k).name,'cell');
            if isempty(id) 
                id = regexp(d(k).name,'.c[0-9].');
                if ~isempty(id)
                    cellnumber = sscanf(d(k).name(id(1)+2:end),'%d');
                end
            end
            if ~isempty(id)
                cid = sscanf(d(k).name(id(1)+4:end),'%d');
                if isempty(filename)
                    filename = d(k).name;
                end
                if ncells == 0 %onlye add these if was not a Cells file
                    cells = [cells cid];
                    fprintf('Checking %s\n',d(k).name);
                    A = LoadExpt([dirpath{did} '/' d(k).name]);
                    H = A.Header;
                    if ~isfield(H,'cellnumber')
                        H.cellnumber = GetCellNumber(A);
                    end
                    if isnan(H.cellnumber) && cellnumber > 0
                        H.cellnumber = cellnumber;
                    end
                    cells = [cells H.cellnumber];
                    D.celltrials{j}(H.cellnumber) = length(A.Trials);
                    D.nspikes{j}(H.cellnumber) = length(cat(1,A.Trials.Spikes));
                    D.durations{j}(H.cellnumber) = expt.value(A,'duration');
                    D.cellfiles{j}{H.cellnumber} = H.loadname;
                    if fitcells
                        res = PlotExpt(A,'frbox','noplot','rcnmin',10,'splisdf');
                        efit = FitExpt(res);
                        if isfield(efit,'RespVar');
                            D.respvar{j}(H.cellnumber) = efit.RespVar;
                        elseif isfield(efit,'di')
                            D.di{j}(H.cellnumber) = efit.di;                            
                        end
                        D.fits{j}{H.cellnumber} = efit;
                        D.plots{j}{H.cellnumber} = res;
                    end

                elseif ~ismember(cid, cells)
                    fprintf('%s has Cell %d Not in %s\n',d(k).name,cid, cellsfile);
                    xcells = [xcells cid];
                end
            end
        end
        cells = unique(cells(cells>0));
        D.cells{j} = cells;
        D.xcells{j} = xcells;
        D.ncells(j) = length(unique(cells));
        D.filename{j} = filename;
        nt = E(e).trialcounts(did,suffid);
        if isempty(nt)
            D.trialsrun(j) = 0;
            nt = 0;
        else
            D.trialsrun(j) = nt;
        end
        if isempty(cells)
            D.rates{j} = 0;
            if D.nblock(j) > 0
                D.trialsrun(j) = nt;
            end
            D.ntrials(j) = 0;
        else
            D.rates{j}(cells) = D.nspikes{j}(cells)./(D.celltrials{j}(cells) .* D.durations{j}(cells));
            D.ntrials(j) = nt;
        end
        if D.nfiles(j) == 0 && D.nblock(j) > 0
            fprintf('%s no %s Expt Files but %d blocks were run (%d Trials)\n',dirpath{did},suffix,D.nblock(j),nt);
            Comments = expt.AddComments(dirpath{did});
            if isfield(Comments,'comment')
                id = find(CellToMat(regexp({Comments.comment},suffix)));
                for k = 1:length(id)
                    cprintf('blue','%s:%s\n',dirpath{did},Comments(id(k)).comment);
                    D.Comments{j}{k} = Comments(id(k)).comment;
                end
            end
        else
            fprintf('%s %d Expt Files (%d cells) match %s,  %d blocks were run (%d Trials)\n',dirpath{did},length(d),D.ncells(j),suffix,D.nblock(j),nt);
        end
    end
    DiskFiles(e) = D;
end

function [names, ids]= FindExpts(namelist, varargin)
names = {};
details = [];

if isempty(namelist)
    ids = [];
    return;
end
exact = 0;

j = 1;
while j <= length(varargin)
    if iscell(varargin{j}) && isfield(varargin{j}{1},'dirpath')
        exdetails = varargin{j};
    elseif strcmp(varargin{j},'exact')
        exact = 'exact';
    end
    j = j+1;
end

if isstruct(namelist{1}) && isfield(namelist{1},'dirpath') %just details
    exdetails = namelist;
    namelist = {};
    for j = 1:length(exdetails)
        if isfield(exdetails{j},'names') %actually has data
            namelist{j} = exdetails{j}.names;
        end
    end
end

if iscellstr(namelist) && isdir(namelist{1}) %cell array of dir names
    for j = 1:length(namelist)
        [names{j}, ids{j}] = ListExptDir(namelist{j},varargin{:});
    end
    return;

elseif iscell(namelist{1}) %not a cell array of strings = array of cellstrings
    ids = {};
    for j = 1:length(namelist)
        if ~isempty(namelist{j})
            good(j) = 1;
        else
            good(j) = 0;
        end
    end
    if cellstrcmp('empty',varargin);
        exdetails = CellToStruct(exdetails(find(good==0)));
        names = {exdetails.dirpath};
        ids = exdetails;
        return;
    end
    namelist = namelist(find(good)); %need to remove empties
    exdetails = exdetails(find(good));
    goodlist = find(good);
    n2 = [];
    matchcounts = [];
    trialcounts = [];
    for j = 1:length(namelist)
%names is a list of all blocks in namelist that match        
        [names{j}, id{j}] = FindExpts(namelist{j},varargin{:},exact);
        details.matches(j) = length(names{j});
        if j <= length(exdetails) && isfield(exdetails{j},'ntrials')
            details.ntrials(j) = sum(exdetails{j}.ntrials(id{j}));
            details.filename{j} = unique(exdetails{j}.filename(id{j}));
        end
        details.matchnames{j} = unique(names{j});
        for k = 1:length(details.matchnames{j})
            [a,eid] = FindExpts(namelist{j},details.matchnames{j}{k},'exact');
            matchcounts(j,k) = length(a);
            trialcounts(j,k) = sum(exdetails{j}.ntrials((eid)));
            if isfield(exdetails{j},'n2')
                n2(j,k) = max(exdetails{j}.n2(eid));
                nt(j,k) = max(exdetails{j}.nt(eid));
                minn2(j,k) = min(exdetails{j}.n2(eid));
                minnt(j,k) = min(exdetails{j}.nt(eid));
            end
        end
    end
    details.ids = id;
    clear id;
    %E = CellToStruct(exdetails);
    id = find(details.matches);
    details.ids = details.ids(id);
    details.matches = details.matches(id);
    details.ntrials = details.ntrials(id);
    details.filename = details.filename(id);
    details.matchnames = details.matchnames(id);
    details.matchid = goodlist(id);
    details.names = names(id);
    details.matchcounts = matchcounts(id,:);
    details.trialcounts = trialcounts(id,:);
    if ~isempty(n2)
    for j = 1:length(id)
        details.n2 = n2(id,:);
        details.nt = nt(id,:);
    end
    end
    allnames = {};
    for j = 1:length(details.matchnames)
        allnames = {allnames{:} details.matchnames{j}{:}};
    end
    [details.namecounts,details.allnames] = Counts(allnames);
    if ~isempty(details)
        for j = 1:length(id)
            details.dirpath{j} = exdetails{id(j)}.dirpath;
        end
    end
    ids = names;
    names = details;
    return;
end


%get here if its a cell array of strings that are not dir names
%= should be a list of expt names
if strcmp(exact,'exact')
    exact = 1;
end

findstr = {};
ids = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'test',4)
    elseif ischar(varargin{j})
        findstr = {findstr{:} varargin{j}};
    end
    j = j+1;
end
nf = 0;
for k = 1:length(findstr)
for j = 1:length(namelist)
    if (~isempty(strfind(namelist{j},findstr{k})) && exact == 0) || ...
            strcmp(namelist{j},findstr{k})
        nf = nf+1;
        names{nf} = namelist{j};
        ids(nf) = j;
    end
end
end


function s = MakeName(E,j)


Expt.Stimvals.et = E.params{j}{1};
if length(E.params{j}) > 1
    Expt.Stimvals.e2 = E.params{j}{2};
else
    Expt.Stimvals.e2 = 'e0';
end
if length(E.params{j}) > 2
    Expt.Stimvals.e3 = E.params{j}{3};
else
    Expt.Stimvals.e3 = 'e0';
end
Expt.Stimvals.st = 'rds';
s = Expt2Name(Expt);
s = strrep(s,'rds.','');
names = {'or' 'orXob' 'sOXorXor' 'orXme'};
suffixes = {'OT' 'ORBW' 'OPP' 'OXM'};
id = find(strcmp(s,names));
if ~isempty(id)
    s = suffixes{id};
end

function [a,b] = CountReps(E)
et = E.Stimvals.et;
e2 = E.Stimvals.e2;
e3 = E.Stimvals.e3;
E.Stimvals.e0 = 0;
E = FillTrials(E,'e0');
x = [E.Trials.(et)];
y = [E.Trials.(e2)];
z = [E.Trials.(e3)];

xv = unique(x);
yv = unique(y);
zv = unique(z);
b = length(xv).*length(yv).*length(zv);
a = length(E.Trials)./b;


function [names, details] = ListExptCells(E, varargin)
showvals = [];
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'show',5)
        j = j+1;
        showvals{end+1} = varargin{j};
    end
    j = j+1;
end

for j = 1:length(E)
    if ~isempty(E{j})
    if ~isfield(E{j}.Header,'expname')
        E{j}.Header.expname = Expt2Name(E{j});
    end
    names{j} = E{j}.Header.expname;
    if ~isfield(E{j}.Stimvals,'ei')
        E{j}.Stimvals.ei = 0;
    end
    if isfield(E{j},'Trials')
        details.ntrials(j) = length(E{j}.Trials);
    elseif isfield(E{j},'ntrials')
        details.ntrials(j) = E{j}.ntrials;
    end
    fn = fieldnames(E{j}.Trials);
    fn = setdiff(fn,{'Start' 'TrialStart' 'End' 'dur' 'id' 'TrueEnd' 'bstimes' 'delay' 'Trial' 'exvals' 'endevent' 'Spikes' 'count' 'OptionCode' 'rwtimes'});
    fn = setdiff(fn,{'ch' 'op' 'rw' 'se' });
    str = [];
    idstr = [];
    if isfield(E{j}.Trials,'id')
        idstr = sprintf(' Id%d-%d',minmax([E{j}.Trials.id]));
    end
    
    for k = 1:length(fn)
        str = [str fn{k} ','];
    end
    xstr = ExptString(E{j},showvals);
    for k = 1:length(showvals)
        details.(showvals{k})(j) = GetEval(E{j},showvals{k});
    end
    if isfield(E{j}.Header,'human')
        [a,b] = CountReps(E{j});
        fprintf('%d: (%d * %.1frep = %dTrials) %s\n',j,b,a,details.ntrials(j),Expt2Name(E{j}));
    else
        fprintf('%d: (%dTrials%s) %s %s %s. Fields %s\n',j,details.ntrials(j),idstr,names{j},xstr,Expt2Name(E{j}),str);
    end
    end
end


function [names, details] = ListExptDir(E, varargin)


relist = 0;
listfullv = 0;
checkupdate = 0;
verbose = 1;
j = 1;
while j <=length(varargin)
    if strncmpi(varargin{j},'check',5)
        checkupdate = 1;
        verbose = 0;
    elseif strncmpi(varargin{j},'fullv',5)
        listfullv = 1;
    elseif strncmpi(varargin{j},'rebuild',5)
        relist = 2; %make list again, but don't call APlaySpkFile
    elseif strncmpi(varargin{j},'relist',5)
        relist = 1;
    elseif strncmpi(varargin{j},'silent',5)
        verbose = 0;
    elseif strncmpi(varargin{j},'verbose',5)
        verbose = 2;
    end
    j = j+1;
end

details.dirpath = E;
if listfullv
    expts = [];
    probes = [];
    sizes = [];
    d = dir([E '/Expt*FullV.mat']);
    for j = 1:length(d)
        expts(j) = GetExptNumber(d(j).name);
        probes(j) = GetProbeFromName(d(j).name);
        sizes(j) = d(j).bytes;
    end
    names = unique(expts);
    details.expts = expts;
    details.probes = probes;
    details.sizes = sizes;
    return;
end

outname = [E '/AllExptList.mat'];
if exist(outname) && relist == 0
    d = dir(outname);
    listdate = d.date;
    load(outname);
    done = 1;
    names = details.names;
    if checkupdate && ~isfield(details,'nt') %old format
        [names, details] = ListExptDir(E, varargin{:},'rebuild');
        return;
    end
    if isfield(details,'listdate')
        listdate = details.listdate;
    end
    if ~isfield(details,'type') && isfield(details,'filename')
        nf = length(unique(details.filename));
        if nf == 1
            details.type = 'WaveMark';
        elseif nf > length(details.filename)/2
            details.type = 'Continuous';
        else
            details.type = 'Unknown';
        end
    end
    if ~isfield(details,'dirpath')
        details.dirpath = E;
    end
    if ~isfield(details,'exptno') || isempty(details.exptno)
        done = 0;
        if verbose > 0
            fprintf('%s no Expts when Last listed\n',E);
        end
    elseif isfield(details,'filename') && checkupdate
        filenames = unique(details.filename);
        for j = 1:length(filenames)
            d = dir([E '/' filenames{j}]);
            if isempty(d)
                fprintf('%s has gone!!\n',filenames{j});
            elseif d.date > listdate
                fprintf('%s has Changed - relisting\n',filenames{j});
                done = 0;
            end
        end
%On Windows the dir call in NOT case sensitive        
        d = dir([E '/*idx.mat']);
        gid = find(CellToMat(regexp({d.name},'.*idx.mat')));
        d = d(gid);
        new = setdiff({d.name},details.filename);
        for j = 1:length(new)
            name = [E '/' new{j}];
            X = load(name);
            if isfield(X,'ExptList') && ~isempty(X.ExptList)
                fprintf('ListExpts: %s is  new\n',new{j});
                done = 0;
            end
        end
    elseif verbose
        for nx = 1:length(names)
            if isfield(details,'goodtrials')
                fprintf('%d:%s %d Trials X %.1f sec\n',details.exptno(nx),names{nx},details.goodtrials(nx),details.trialdur(nx)./10000);
            else
                fprintf('%d:%s %d Trials\n',details.exptno(nx),names{nx},details.ntrials(nx));
            end
        end
    end
    if done
        return;
    end
end

if verbose > 0
fprintf('Listing %s\n',E);
end
d = dir([E '/*idx.mat']);
if ~isempty(d)
    exid = GetExptNumber({d.name});
    [~, sortid] = sort(exid);
    d = d(sortid);
end
nx = 0;
names = {};
details = [];
details.dirpath = E;
for j = 1:length(d)
    if sum(strcmp(d(j).name,{'FileIdx.mat' 'OldIdx.mat'}))
    else
    name = [E '/' d(j).name];
    clear ExptList;
    load(name);
    id = regexp(d(j).name,'\.[0-9]*idx');
    if isempty(id)
        fileno = 0;
    else
        fileno = sscanf(d(j).name(id(1)+1:end),'%d');
    end
    if ~exist('ExptList','var') && exist('Expt','var')
        ExptList = BuildExptList(Expt, Expts);
    end
    if exist('ExptList','var') && isfield(Expts,'result')
        for k = 1:length(ExptList)
            Ex = ExptList(k);
            e = find([Expts.start] == ExptList(k).start);
            if isempty(Expts(e).result)
                Expts(e).result = 2;
            end
            if length(e) && Expts(e).result == 2
                nx = nx+1;
                names{nx} = ExptList(k).expname;
                trials = Expts(e).firsttrial:Expts(e).lasttrial;
                trials(Expt.Trials.Result(trials)==0) = [];
                trials(trials > length(Expt.Trials.Start)) = [];
                tids = Expt.Trials.id(trials);
                details.trialdur(nx) = mean(Expt.Trials.End(trials) - Expt.Trials.Start(trials));
                details.ntrials(nx) = Expts(e).lasttrial-Expts(e).firsttrial;
                details.goodtrials(nx) = length(trials);
                details.filename{nx} = d(j).name;
                details.exptno(nx) = k+fileno-1;
                if isempty(tids)
                    details.idrange(nx,:) = [0 0];
                else
                    details.idrange(nx,:) = minmax(tids);
                end
                exs = {};
                exs{1} = Ex.et;
%need to update this to use exptvars for manual expts one day...                
                if isfield(Ex,'e2') && ~strcmp(Ex.e2,'e0')
                    exs = {exs{:} Ex.e2};
                end
                if isfield(Ex,'e3') && ~strcmp(Ex.e3,'e0')
                    exs  = {exs{:} Ex.e2};
                end
                details.params{nx} = exs;
                if isfield(Expt.Trials,exs{1})
                    if iscell(Expt.Trials.(exs{1}))
                        [a,b] = Counts(cat(1,Expt.Trials.(exs{1}){:}));
                    else
                        [a,b] = Counts(Expt.Trials.(exs{1}));
                    end
                elseif strcmp('sz',exs{1})
                    [a,b] = Counts(Expt.Trials.wi);
                else
                    a = 1;
                end
                details.nt(nx) = sum(a>=mean(a)/2);
%Check for missng fields we need, but not if they can be Filled                    

                if sum(strcmp('ce',exs)) && isfield(Expt.Trials,'ce')
                    if iscell(Expt.Trials.ce)
                        [a,b] = Counts(cat(1,Expt.Trials.ce{:}));
                    else
                        [a,b] = Counts(Expt.Trials.ce);
                    end
                    details.n2(nx) = sum(a>mean(a)/2);
                elseif sum(strcmp(Ex.e2,{'e0' ''}))
                        details.n2(nx) = 0;
                elseif ~isfield(Expt.Trials,Ex.e2) 
                    fillf = {'cL' 'cR'};
                    if sum(strcmp(Ex.e2,fillf)) ==0
                        details = AddError(details,'-show','ListExpts: No Field for Expt %s',Ex.e2);
                        details.n2(nx) = 0;
                    else
                        id = find(strcmp(Ex.e2,fillf));
                        if ismember(id,[1 2]) && isfield(Expt.Trials,'ic')
                            [a,b] = Counts(Expt.Trials.ic);
                            details.n2(nx) = sum(a>mean(a)/2);
                        else
                        details.n2(nx) = 0;  %might want to make this fill and count later..
                        end
                    end
                else
                    [a,b] = Counts(Expt.Trials.(Ex.e2));
                    details.n2(nx) = sum(a>mean(a)/2);
                end
                if details.nt(nx) > 20 || details.n2(nx) > 20
                    fprintf('*');
                end
                if verbose
                fprintf('%d:%s %d Trials X %.1f sec %d X %d stim\n',details.exptno(nx),names{nx},details.goodtrials(nx),details.trialdur(nx)./10000,details.nt(nx),details.n2(nx));
                end
            end
        end
    end
    end
end
if nx > 0
    [a,b] = sort(details.exptno);
    names = names(b);
    if sum(fileno) == 0
        details.type = 'WaveMark';
    elseif sum(fileno > 0) > length(a)/2
        details.type = 'Continuous';
    else
        details.type = 'Unknown';
    end
    details.ntrials = details.ntrials(b);
    details.filename = details.filename(b);
    details.exptno = details.exptno(b);
    details.names = names;
    try
        save(outname,'details');
    end
elseif ~exist(outname) || relist == 1
    if ~isempty(strfind(E,'icarus'))
             [a,b] = fileparts(E);
        mfile = [E '/ic' b '.mat'];
    else
        [a,b,c,d] = GetMonkeyName(E);
        mfile = [E '/'  b c '.mat'];
    end
    details.ntrials = [];
    details.filename = {};
    details.exptno = [];
    details.type = 'Empty';
    details.names = {};
    if exist(mfile) && relist ~= 2
        fprintf('Building idx for %s\n',mfile);
        APlaySpkFile(mfile);
        [names, details] = ListExptDir(E,'rebuild');
        if ~exist(outname)
            save(outname,'details');
        end
        return;
    end
    save(outname,'details');
else
    details.ntrials = [];
    details.filename = {};
    details.exptno = [];
    details.type = 'Empty';
    details.names = {};
end


function ExptList = BuildExptList(Expt, Expts)

ExptList = [];

nx = 0;
for j = 1:length(Expts)
    if Expts(j).lasttrial - Expts(j).firsttrial > 10
        nx = nx+1;
        ts = Expts(j).firsttrial;
        if isfield(Expt.Trials,'et')
            ExptList(nx).et = Expt.Trials.et{ts};
            ExptList(nx).e2 = Expt.Trials.e2{ts};
            ExptList(nx).e3 = Expt.Trials.e3{ts};
        elseif isfield(Expts,'et')
            ExptList(nx).et = Expts(j).et;
            ExptList(nx).e2 = Expts(j).e2;
            ExptList(nx).e3 = Expts(j).e3;
        end
        ExptList(nx).start = Expts(j).start;
        ExptList(nx).expname = [ExptList(nx).et 'X' ExptList(nx).e2 'X' ExptList(nx).e3];
        if isfield(Expt.Trials,'OptionCode')
        if strfind(Expt.Trials.OptionCode{ts},'+fS')
            ExptList(nx).expname = [ExptList(nx).expname 'FS'];
        end
        end
    end
end

function [names, details] = ListPsychDir(E, depth, varargin)
names = {};
details = {};
d = dir(E);
if isdir(E)
    root = [E '/'];
else
    id = strfind(E,'/');
    if length(id)
        root = E(1:id(end));
    end
end

nd = 0;
ngood = 0;

for j = 1:length(d)
     path = [root d(j).name];
     [a,b,c] = fileparts(d(j).name);
     if ~d(j).isdir && isempty(c)
     expts = PsychMon(path,'getexpts');
     if ~isempty(expts)
     [names{j}, details{j}] = ListExpts(expts);
     end
     end
end


function [names, details] = ListExptDirs(E, depth, varargin)

names = {};
details = {};
d = dir(E);
if isdir(E)
    root = [E '/'];
else
    id = strfind(E,'/');
    if length(id)
        root = E(1:id(end));
    end
end

nd = 0;
ngood = 0;

for j = 1:length(d)
    if d(j).isdir && d(j).name(1) ~= '.'
        if depth > 1
            path = [root d(j).name];
            [a, b] = ListExptDirs(path, depth-1, varargin{:});
            for k = 1:length(a)
                nd = nd+1;
                names{nd} = a{k};
                details{nd} = b{k};
                ngood = ngood+1;
            end
        else
            nd = nd+1;
            path = [root d(j).name];
            [names{nd}, details{nd}] = ListExptDir(path, varargin{:});
            details{nd}.dirpath = path;
            if ~isempty(names{nd})
                ngood = ngood+1;
            end
        end
    end
end

if ngood == 0
    names = {};
    details = {};
end