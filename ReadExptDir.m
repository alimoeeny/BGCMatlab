function [AllExpts, details, Ex] = ReadExptDir(name, varargin)
% [Expts, details, Idx] = ReadExptDir(name, varargin)
% Read all SPike2 .mat files in a directory and combines the Expts lists
% into one list. Saves the result. nameExpts.mat
% ReadExptDir(name, 'relist') Calls APlaySpkFile and rebuilds all the
% idx.mat files
% ReadExptDir(name, 'resort') rebuilds Expt.Trials from the idx file
%        Use this also to check for new .mat files that have not yet had
%        APlaySpkFile run.
%
% ReadExptDir(name, 'fixpen') rebuilds the .ufl files recording penetration
%                             number and RF position
%       These both call SetTrialOffsets, which forces monotonically
%       increasing Expt.Trials.Trial across all blocks
% 
% ReadExptDir(name, 'showerrs') prints out errors warnings from
% APlaySpkFile/CheckExpts
%
% ReadExptDir(name, 'checkstimtime') Checks duration in serial output file
% Called by AplaySpkFile  if first argument is a directory
%Idx returns the contents of the Idx files made by APlaySpkFile. 
%      requesting this adds considerably to the load time for previously
%      indexed folders
% 
% ReadExptDir(name, 'addmains') adds timestamps for frame pulse and mains
% phase
%
%                ...,'summary') Just loads a summary file without
%                individual trials - gets a quick list of expts...

state.relist = 0;
state.resort = 0; %redo SortExpts, but not whole listing
state.online = 0;
state.quick = 0;
state.summary = 0;
state.checkstimtime =0;
state.checkdur=0;
state.checkframe=0;
state.showerrs = 0;
state.mkufl = 0;
state.addmains = 0;
state.verbose = [0 1 1];


j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'addmains',4)
        state.addmains=1;
    elseif strncmpi(varargin{j},'online',4)
        state.online=1;
    elseif strncmpi(varargin{j},'resort',4)
        state.resort=1;
    elseif strncmpi(varargin{j},'checkall',4)
        state.checkdur=1;
        state.checkframe=1;
    elseif strncmpi(varargin{j},'fixpen',4)
        state.mkufl = 1;
    elseif strncmpi(varargin{j},'checkstimstime',4)
        state.checkstimtime=1;
    elseif strncmpi(varargin{j},'quiet',4)
        state.verbose = [0 0 1];
    elseif strncmpi(varargin{j},'silent',4)
        state.verbose = [0 0 0];
    elseif strncmpi(varargin{j},'silent',4)
        state.verbose = [1 1 1];
    elseif strncmpi(varargin{j},'quick',4)
        state.quick=1;
    elseif strncmpi(varargin{j},'summary',4)
        state.summary=1;
    elseif strncmpi(varargin{j},'showerrs',6)
        state.showerrs =1;
    elseif strncmpi(varargin{j},'relist',4)
        state.relist = 1;
        state.resort = 1;
    end
    j=j+1;
end


if ~isdir(name) && isdir(fileparts(name))
    filename = name;
    name = fileparts(filename);
end

expdates = [];
%first sort numerically by suffix number
%And find modification times of individual expts files
d = mydir([name '/*.mat']);
if isempty(d)
    AllExpts = {};
    Ex = {};
    cprintf('red','No Matlab files in %s\n',name);
    return;
end
mnk =GetMonkeyName(name);
suffixes = [];
expfiles = {};
matfiles = {};
AllErrs = {};
BinocFile = '';
details.needfullv = [];

if state.summary 
    outname  = BuildFileName(name,'ExptSummary');
    if exist(outname,'file')
        load(outname);
        AllExpts = Expts;
        return;
    else
    end
end

for j = 1:length(d)
    if state.online
        if ~isempty(regexp(d(j).filename,'Expt[0-9]*.mat'))
            suffixes(j) = str2double(regexprep(d(j).filename,'Expt([0-9]*).mat','$1'));
            truesuffix(j) =1;
        end
    elseif ~isempty(regexp(d(j).name,[mnk '[E,M,\.,G,0-9]*\.[0-9]+.mat']))
        suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z,0-9]\.([0-9]+).mat','$1'));        
        truesuffix(j) = 1;
        BinocFile = regexprep(d(j).name,'\.[0-9]*.mat','.bnc.mat');
        matfiles{suffixes(j)} = d(j).name;
    elseif ~isempty(regexp(d(j).name,[mnk '[E,M,\.,G,0-9]+.mat'])) %not a .N.mat file, but a XXX.mat file 
%get suffix from first part of file name. Beware of extra files in a folter
%like jbeG068
        [a,b,c] = GetMonkeyName(d(j).name);
        if strfind(d(j).filename,c)
            suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z]([0-9]+).mat','$1'));
        else
            suffixes(j) = 3000 + str2double(regexprep(d(j).filename,'.*[\.,A-z]([0-9]+).mat','$1'));
        end
        truesuffix(j) = 0;
    elseif ~isempty(regexp(d(j).name,[mnk '[E,M,\.,G,0-9]*Expts.mat']))
        e = GetExptNumber(d(j).name);
        if e > 0
            expdates(e) = d(j).datenum;
            expfiles{e} = d(j).name;
        end
    elseif ~isempty(regexp(d(j).name,[mnk 'Expt[0-9]+FullV.mat']))
        e = GetExptNumber(d(j).name);
        needfullv(e) = 0;
    end
end

if isempty(suffixes) %if no suffixes, must be in single files. 
    ne = 0;
    for j = 1:length(d)
        if ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*idx.mat']))
            e = GetExptNumber(d(j).name);
            ne = ne+1;
            suffixes(j) = ne;
        end
    end
end

[a,b] =sort(suffixes);
sid = b(a> 0 & a < 2000); %get rid of Utah ns5->.mat files 

if length(sid) == 1 && truesuffix(sid) == 0
    suffixes(sid) = 1;
%    sid = [];
end

[a,b,c] = GetMonkeyName(name);
outname  = [name '/' a c 'Expts.mat'];
errfile = [name '/ClusterErrors.mat'];

if state.relist && exist(BinocFile);
    gotbnc = 0;
    if isappdata(0,'BinocExpt')
        BinocExpt = getappdata(0,'BinocExpt');
        if isfield(BinocExpt,'loadname') && strcmp(BinocExpt.loadname,BinocFile)
            gotbnc = 1;
        end
    end
    if gotbnc == 0
        fprintf('Loading %s\n',BinocFile);
        load(BinocFile);
%    BinocExpt.readdate = now;
         BinocExpt.loaddate = now;
         setappdata(0,'BinocExpt',BinocExpt);
    end
elseif isappdata(0,'BinocExpt')
    rmappdata(0,'BinocExpt');
end


if state.checkstimtime
    CheckStimTimes({d(sid).name});
end
if state.mkufl %rebuild .ufl files
    for j = sid(:)'
        APlaySpkFile(d(j).name,'rfs');
    end
    BuildRFData(name,'rebuild');
    state.resort = 1;
end

if nargout > 2 || state.addmains
    needidx = 1;
else
    needidx = 0;
end
    
if exist(outname) && state.resort == 0
        cerrexlist = [];
        xd = dir(outname);
        ts = now;
        
    load(outname);
    if mytoc(ts) > 5
        fprintf('Loading %s took %.2f = %.1f Mb/sec\n',outname,mytoc(ts),xd.bytes./(mytoc(ts).*1024.*1024));
    end
        
    idxname = strrep(outname,'Expts','ExptIdx');
    if needidx && ~exist('Idx','var') && exist(idxname)
        load(idxname);        
    elseif exist('AllErrs','var')
        details.AllErrs = AllErrs;
    else
        details.AllErrs = {};
    end
    if exist('Idx','var')
        details.gotidx = zeros(size(Expts));
        details.gotidx(1:length(Idx)) = 1;
    end
    if exist(errfile)
        load(errfile);
    else
        errorlist = [];
    end
    if isfield(errorlist,'ex')
        cerrexlist = [errorlist.ex];
    end
    if ~exist('Idx','var') && ~exist(idxname) %old file
        state.resort = 1;
    else
    reloaded = 0;
    for j = 1:length(Expts)
        if ~isfield(Expts{j}.Header,'suffix')
            Expts{j}.Header.suffix = GetExptNumber(Expts{j});
        end
        e = Expts{j}.Header.suffix;
        if e > 0
            errid = find(ismember(cerrexlist,e));
            if ~isempty(errid)
                Expts{j}.clustererrs = errorlist(errid);
            end
            if e > 0 && length(expdates) >= e && expdates(e) > xd.datenum
                fprintf('Reloading %d\n',e);
                E = load(expfiles{e});
                if ~isempty(E.Expts)
                    Expts{j} = E.Expts{1};
                end
                Idx{j} = E.Tidx;
                reloaded = reloaded+1;
            end
            exptgood(e) = e;
        else
            exptgood(j) = j;
        end
        if isfield(Expts{j}.Header,'loadname')
            a= fileparts(Expts{j}.Header.loadname);
           loadname = strrep(Expts{j}.Header.loadname,a,name);
           Expts{j}.Header.loadname = loadname;
        end
        if ~isfield(Expts{j},'errs') && ~isfield(details,'AllErrs') && length(Idx) == length(Expts) && isfield(Idx{j},'errs')
            if ~isempty(Idx{j}.errs)
                Expts{j}.errs = Idx{j}.errs;
                if isfield(Idx{j},'errdata')
                    Expts{j}.errdata = Idx{j}.errdata;
                end
            end
        end
        if state.addmains && isfield(Expts{j},'Trials')
            Expts{j}.mainstimes = Idx{j}.mainstimes;
            Expts{j}.frametimes = Idx{j}.frametimes;
        end
    end
    id = find(exptgood ==0);
    for j = 1:length(id) %Expts not saved in Expts for some reason
        e = id(j);
        if length(expdates) >= e && expdates(e) > 0 
            X = load(expfiles{e});
            if ~isempty(X.Expts)
                cprintf('red','Adding Expt %d (missing) from %s\n',e,expfiles{e});
                Expts(e+1:end+1) = Expts(e:end);
                Expts{e} = X.Expts{1};
                reloaded = reloaded+1;
            else
                cprintf('blue','Expt %d (missing) Looks to be empty. Consider Deleting these files\n',e);
            end
        end
    end
    id = setdiff(suffixes(sid),exptgood);
    for j = 1:length(id)
        if id(j) > length(matfiles)
            fprintf('No .mat file for Suffix %d\n',id(j));
        else
            fprintf('%s Not processed\n',matfiles{id(j)});
        end
    end
    if reloaded
        [a,b, Expts, AllErrs] = CheckExpts(Expts,'quiet');
        fprintf('Saving Expts with reloaded suffixes\n');
        save(outname,'Expts');
        idxname = strrep(outname,'Expts','ExptIdx');
        save(idxname, 'Idx','AllErrs');
    elseif state.checkframe
        [a,b, Expts] = CheckExpts(Expts,'quiet');
    else
        [a,b, Expts] = CheckExpts(Expts,'quick');
    end
    if nargout > 1
        x = CheckFiles(name, Expts);
        details = CopyFields(details,x);
    end
    if isempty(AllErrs)
        AllErrs = {};
        errdata = [];
        for j = 1:length(Expts)
            ne = length(AllErrs);
            eid = GetExptNumber(Expts{j});
            if isfield(Expts{j},'errs')
                if iscellstr(Expts{j}.errs)
                    AllErrs = {AllErrs{:} Expts{j}.errs{:}};
                    if isfield(Expts{j},'errdata')
                        errdata = CatStruct(errdata, Expts{j}.errdata);
                    else
                        [errdata(ne+1:length(AllErrs)).time] = deal(0);
                    end
                    [errdata(ne+1:length(AllErrs)).exptno] = deal(eid);
                elseif isfield(Expts{j}.errs,'errmsg')
                    for k = 1:length(Expts{j}.errs)
                        if iscell(Expts{j}.errs(k).errmsg)
                            AllErrs{end+1} = Expts{j}.errs(k).errmsg{1};
                        elseif ischar(Expts{j}.errs(k).errmsg)
                            AllErrs{end+1} = Expts{j}.errs(k).errmsg;
                        end
                        errdata = CopySFields(errdata,length(AllErrs),Expts{j}.errs(k));
                        errdata(length(AllErrs)).expto = eid;
                    end
                end
            end
        end
        if length(AllErrs) == length(errdata)
            for j = 1:length(errdata)
                errdata(j).s = AllErrs{j};
            end
            details.Errors = errdata;
        else
            details.errs = AllErrs;
            details.errdata = errdata;
        end
        if state.verbose(2)
            if length(errdata) == length(AllErrs)
                ShowErrors(errdata,'severe');
            else
                ShowErrors(AllErrs,'severe');
            end
        end
    end
    details.AllErrs = AllErrs;
    details = FixDetails(details);
    
    AllExpts = SetTrialOffsets(Expts);
    AllExpts = expt.SetTimeOffsets(AllExpts);
    AllExpts = expt.AddComments(AllExpts);
    if nargout > 2
        Ex = Idx;
    end
    if state.showerrs
        details.errdata = ShowExptErrs(AllExpts);
    end
    ListExpts(name,'check');
    if state.summary %get here if 'summary' argument given, but no file
        SaveExptSummary(name, AllExpts);
    end
    return;
    end
end

AllExpts = {};
ExtraExpts = {};
AllErrs = {};
combineexpts = 0;

if isempty(sid) %No files with Expt data
    [a,b] = fileparts(name);
    if strcmp(b,'Ecker')
        [AllExpts, details, Ex] = ReadExptDir(a, varargin{:});
    end
    return;
end
nex = 1;
Default.toplevel = 0;  %let aplayspkfile store app data
for j = 1:length(sid)
    
    fprintf('Reading %s\n',d(sid(j)).name);
    [Ex{nex}, Expts] = APlaySpkFile(d(sid(j)).name,'nospikes','noerrs', 'Defaults', Default, varargin{:});
    if length(Expts) > 1 && length(sid) > 2 %Usually this is an error, and the first expt is bad.  
%Could add a test later to try and combine these if possible        
        fprintf('%d Expts in %s\n',length(Expts),d(sid(j)).name);
        if combineexpts
            AllExpts = {AllExpts{:} CombineExpts(Expts)};
        else
            AllExpts = {AllExpts{:} Expts{end}};
            ExtraExpts = {ExtraExpts{:} Expts{1:end-1}};
        end
    elseif ~isempty(Expts)
        AllExpts = {AllExpts{:} Expts{:}};
    end
    if isfield(Ex{nex},'errs') && ~isempty(Ex{nex}.errs)
        cprintf('blue','Errors for %s\n',d(sid(j)).name);
        for k = 1:length(Ex{nex}.errs)
            cprintf('red','%s\n',Ex{nex}.errs{k});
        end
        if isempty(Expts)
            AllErrs = {AllErrs{:} Ex{nex}.errs};
        elseif isfield(AllExpts{end},'errs')
            if isfield(AllExpts{end}.errs,'msg')
                for k = 1:length(Ex{nex}.errs)
                    AllExpts{end}.errs(end+1).msg = Ex{nex}.errs{k};
                    AllExpts{end}.errs(end+1).t = 0;
                end
            elseif isfield(AllExpts{end}.errs,'errmsg')
                for k = 1:length(Ex{nex}.errs)
                    AllExpts{end}.errs(end+1).errmsg = Ex{nex}.errs{k};
                    AllExpts{end}.errs(end+1).t = 0;
                end
            else
                AllExpts{end}.errs = {AllExpts{end}.errs{:} Ex{nex}.errs{:}};
            end
        else
            AllExpts{end}.errs = Ex{nex}.errs;
        end
    end
    nex = length(AllExpts)+1; %so that it lines up with Expts,
end
AllExpts = SetTrialOffsets(AllExpts);
Array = GetArrayConfig(name);


for j = 1:length(AllExpts)
    e = GetExptNumber(AllExpts{j});
    if isfield(AllExpts{j},'Header') && e > 1
        AllExpts{j}.Header.exptno = e;
        AllExpts{j}.Header.arrayname = str(Array);
    end
    if e > 0 && e <= length(expfiles) && exist(expfiles{e})
        E = load(expfiles{e});
        if isfield(E,'Expts') && ~isempty(E.Expts)
            exptno = GetExptNumber(E.Expts{1});
            if exptno > 1
                E.Expts{1}.Header.exptno = exptno;
            end
            E.Expts{1}.Header.trialoffset = AllExpts{j}.Header.trialoffset;
            save(expfiles{e},'-struct','E');
        end       
    end
end

for j = 1:length(ExtraExpts)
    AllExpts{end+1} = ExtraExpts{j};
    AllExpts{end}.Header.exptno = 100 + GetExptNumber(ExtraExpts{j});
    Ex{length(AllExpts)}.Header = ExtraExpts{j}.Header;
end

for j = 1:length(Ex)
    if isfield(Ex{j},'errs') && ~isempty(Ex{j}.errs)
        cprintf('blue','Errors for %s\n',d(sid(j)).name);
        for k = 1:length(Ex{j}.errs)
            cprintf('red','%s\n',Ex{j}.errs{k});
        end
    end
end
AllExpts = expt.SetTimeOffsets(AllExpts);
AllExpts = expt.fix(AllExpts,'ed');
[a,b,c] = GetMonkeyName(name);
outname  = [name '/' a c 'Expts.mat'];
Expts = AllExpts;
Idx = Ex;

if ~exist(fileparts(outname))
    cprintf('error','Cannot save %s. Folder does not exist\n',outname);
elseif state.quick == 0 || ~exist(outname) %don't write out Expts if didn't load everything
    [a,b, Expts] = CheckExpts(Expts,'quiet');
    save(outname, 'Expts','AllErrs');
    idxname = strrep(outname,'Expts','ExptIdx');
    save(idxname, 'Idx','AllErrs');
end
details = CheckFiles(name, AllExpts);

if nargout > 1 %? combine Ex{}
    
end

function details= CheckFiles(dirname, Expts)

details.needfullv = [];
for j = 1:length(Expts)
    e = GetExptNumber(Expts{j});
    name = sprintf('%s/Expt%dFullV.mat',dirname,e);
    if exist(name)
        details.needfullv(e) = 0;
    else
        details.needfullv(e) = 1;
    end
    name = sprintf('%s/Expt%dClusterTimes.mat',dirname,e);
    if exist(name)
        details.needclusters(e) = 0;
    else
        details.needclusters(e) = 1;
    end
    name = sprintf('%s/Expt%dAutoClusterTimes.mat',dirname,e);
    if exist(name)
        details.needautoclusters(e) = 0;
    else
        details.needautoclusters(e) = 1;
    end
end

function SaveExptSummary(name, E)

outname = sprintf('%s/ExptSummary.mat',name);
if exist(outname) && 0
    return;
end
try
Stimvals = [];
Stimdifs = {};
for j = length(E):-1:1
    if isfield(E{j},'Stimvals') 
        if isempty(Stimvals)
            Stimvals = E{j}.Stimvals;
            Stimdifs{j} = E{j}.Stimvals;
        else
            f = fields(E{j}.Stimvals);
            for k = 1:length(f)
                if ~isfield(Stimvals,f{k})
                    Stimvals.f{k} = E{j}.Stimvals.(f{k});
                elseif ischar(E{j}.Stimvals.(f{k})) && ~strcmp(E{j}.Stimvals.(f{k}),Stimvals.(f{k}));
                    Stimdifs{j}.(f{k}) = E{j}.Stimvals.(f{k});
                elseif E{j}.Stimvals.(f{k}) ~= Stimvals.(f{k})
                    Stimdifs{j}.(f{k}) = E{j}.Stimvals.(f{k});
                end
            end
        end
    end
end
for j = 1:length(E)
    Expts{j} = rmfields(E{j},{'Trials' 'Stimvals' 'errs'});
    if j <= length(Stimdifs) && ~isempty(Stimdifs{j})
        Expts{j}.Stimvals = Stimdifs{j};
    end
end
%Expts = CellToStruct(Expts);
save(outname,'Expts');
catch ME
    fprintf('Error Trying to save summary file %s\n',outname);
    CheckExceptions(ME);
end
function CheckStimTimes(name)

if iscellstr(name)
    for j = 1:length(name)
        CheckStimTimes(name{j});
    end
    return;
end

aname = regexprep(name,'\.[0-9]+\.mat','A$0');
X = load(name);
f = fields(X);
for j = 1:length(f)
    V = X.(f{j});
    if isfield(V,'title') && strcmp(V.title,'StimOn')
        stimtimes = V.times;
    end
end
A = load(aname);
f = fields(A);
for j = 1:length(f)
    V = A.(f{j});
    if isfield(V,'title') && strcmp(V.title,'StimOn')
        astimtimes = V.times;
    end
end
lendiff = length(astimtimes) - length(stimtimes);
td = NaN;
if lendiff == 0
    xc = corrcoef(astimtimes, stimtimes);
    if xc(1,2) > 0.99
        td = mean(astimtimes-stimtimes);
    end
elseif lendiff == 1
    xc = corrcoef(astimtimes(2:end), stimtimes);
    if xc(1,2) > 0.99
        td = mean(astimtimes(2:end)-stimtimes);
    end
end
if abs(td) > 0.01
    fprintf('TIme mismatch %.3f in %s\n',td, name);
end
function D = FixDetails(D)

if isfield(D,'AllErrs') && ~isempty(D.AllErrs) && ~iscellstr(D.AllErrs)
    errdata = [];
    if isfield(D.AllErrs,'errs')
        for j = 1:length(D.AllErrs.errs)
            if ischar(D.AllErrs.errs{j})
                A{j} = D.AllErrs.errs{j};
            elseif iscellstr(D.AllErrs.errs{j})
                A{j} = D.AllErrs.errs{j}{1};
            else
                A{j} = '';
            end
        end
        D.AllErrs = A;
    end
    if isfield(D,'errdata') && length(D.errdata) == length(D.AllErrs)
        for j = 1:length(D.errdata)
            D.errdata(j).s = D.AllErrs{j};
        end
    end
end