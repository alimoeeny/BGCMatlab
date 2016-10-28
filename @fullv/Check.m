function [result, FullV] = Check(FullV, varargin)
% fullv.Check(FullV, varargin)
% fullv.Check(FullV, Expt)
%Check Integrity of data in FullV Files. Looks for epochs were there are
%missing channels, based on low amplitudes in probes1-16, or 17-24% CheckFullV(FullV, varargin)
%Check Integrity of data in FullV Files. Looks for epochs were there are
%missing channels, based on low amplitudes in probes1-16, or 17-24
%
%If FullV is a directory name, checks all FullVs in the folder
%if FullV is a cell string array, checks each file/folder in arrar
%If FullV is a cell array of FullVs checks each
%
%fullv.Check(dirname,'clearspk')  removes spkblk files IF no errors were detected
%
%fullv.Check(dirname,'clearerr') removes old errors from the FullV files, 
% and removes THESE errors  from errors.mat. So if checking the actual
% FullV files in OK, this should be safe. 
%
%only do this if you are sure old errors are fixed.
%                  ...,'listonly') lists fullV errors in Errors.mat (no new checks)
%                  ...,'listall') lists all errors in Errors.mat
%
%Any Errors are recorded in Errors.mat 

%Old CheckFullV now calls fullv.Check
%"Low" amplitudes are are calculated with SD AFTER adding back in the MeanV
%also check for suspicicous high correlations with meanV

ispk = 22;
scales = 1;
spkrate = 100;
showth = 0;
filtershape = 0;
maxsize = NaN;
fix.scale = 0;
fix.save = 0;
fix.int = 0;
Expt = [];
exptlist = [];
Expts = {};
%DO NOT add fields to result here - breaks references to result{j} =
%fullv.check below, when processing multiple files
clearspk = 0;
clearerrs = 0;
listonly = 0;
forceread = 0;
EX.excluded = {};
checkmains = 0;
verifyerrors = 0; 

sargs = {};
j = 1;

while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'Header')
            Expt = varargin{j};
        elseif isfield(varargin{j},'excluded')
            EX = varargin{j};
        end
    elseif iscell(varargin{j})
        if nansum(CellToMat(varargin{j},'Header.frameperiod'))
            Expts = varargin{j};
        end
    elseif strncmpi(varargin{j},'clearspk',7) %remove spkblk files if no errors
        clearspk = 1;
    elseif strncmpi(varargin{j},'clearerr',7) %remove spkblk files if no errors
        clearerrs = 1;
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        exptlist = varargin{j};
    elseif strncmpi(varargin{j},'fixint',6)
        fix.int =1;
    elseif strncmpi(varargin{j},'fixscale',7)
        fix.scale =1;
    elseif strncmpi(varargin{j},'force',5)
        forceread = 1;
    elseif strncmpi(varargin{j},'listall',7)
        listonly = 2; %just read error files
    elseif strncmpi(varargin{j},'listmatch',7)
        listonly = 3; %just read error files
        j = j+1;
        matcherr = varargin{j};
    elseif strncmpi(varargin{j},'listonly',7)
        listonly = 4; %just read error files
    elseif strncmpi(varargin{j},'list',4)
        listonly = 1; %just read error files
    elseif strncmpi(varargin{j},'mains',5)
        checkmains = 1;
    elseif strncmpi(varargin{j},'maxsize',5)
        j = j+1;
        maxsize = varargin{j};
    elseif strncmpi(varargin{j},'odd',3)
        filtershape = 1;
    elseif strncmpi(varargin{j},'scale',4)
        j = j+1;
        scales = varargin{j};;
    elseif strncmpi(varargin{j},'tchan',4)
        j = j+1;
        ispk = varargin{j};
    elseif strncmpi(varargin{j},'verify',4)
        verifyerrors = 1;
    else
        sargs = {sargs{:} varargin{j}};
    end
    j = j+1;
end


if iscellstr(FullV)
    result = {};
    names = FullV;
    if ~isempty(Expts)
        expts = GetExptNumber(Expts);
    else
        expts = [];
    end
    loadtime = [];
    for j = 1:length(names)
        if isdir(names{j})
            result{j} = fullv.Check(names{j},varargin{:});
        else
            go = 1;
            if maxsize > 0
                d = dir(names{j});
                if d.bytes > maxsize
                    fprintf('Ignoring %s size is %d\n',names{j},d.bytes);
                    go = 0;
                end
            end
            e = GetExptNumber(names{j});
            if e > 0
                if ~isempty(exptlist) && ~ismember(e,exptlist)
                    xstr = '';
                    if e <= length(EX.excluded) && ~isempty(EX.excluded{e})
                        xstr = [' (' EX.excluded{e} ')'];
                    end
                    fprintf('Not Checking %s%s\n',names{j},xstr);
                    go = 0;
                elseif isfield(EX,'needcheck') && ~ismember(e,EX.needcheck)
                    if length(EX.excluded) >= e
                        fprintf('Not Checking %s (%s)\n',names{j},EX.excluded{e});
                    else
                        fprintf('Not Checking %s\n',names{j});
                    end
                    go = 0;
                end
            else
                fprintf('Not Checking %s Cannot determine Expt Number\n',names{j});
                go = 0;
            end
            if go
                a = find(expts == e);
                if ~isempty(a)
                    Expt = Expts{a};
                else
                    Expt = [];
                end
                if exist(names{j})
                    fprintf('Loading %s at %s\n',names{j},datestr(now));
                    FullV = LoadFullV(names{j});
                    loadtime(j) = FullV.loadtime;
                    loadsize(j) = FullV.size;
                    result{j} = fullv.Check(FullV,Expt,varargin{:});
                    result{j} = AddError(result{j},FullV,'-prog','LoadFullV');
                    clear FullV;
                else
                    result{j}.progname = 'fullv.Check';
                    result{j}.eid = e;
                    result{j} = AddError(result{j},'Cannot Read %s: no file',names{j});
                end
            end
        end
    end
    if ~isempty(loadtime)
        readrate = sum(loadsize)./sum(loadtime);
        fprintf('Read Rate %.3f Mb/sec at %s\n',readrate,datestr(now));
    end
    return;
elseif iscell(FullV) %array of results
    if ~isempty(CellToMat(FullV,'checktimes'))
        for j = 1:length(FullV)
            result{j} = ShowResults(FullV{j},varargin{:});
        end
        return;
    elseif sum(CellToMat(FullV,'checktime'))
        ShowResults(FullV,varargin{:});
        return;
    end
elseif ischar(FullV)
    if isdir(FullV)
        ExptErrs = [];
        errfile = [FullV '/Errors.mat'];        
        if exist(errfile)
            X = my.load(errfile);
            if iscell(X.Errors)
                X.Errors = CellToStruct(X.Errors);
            end
            result = X;
            result.name = X.myload.name;
        else
            X.Errors = {};
            result.Errors = [];
        end
        if ~isfield(result,'Verrs')
            result.Verrs = [];
        end
        if listonly
            if exist(errfile)
                if listonly == 2
                   [result.Errors, result.Verrs] = ShowResults(result,'-all');
                elseif listonly == 3
                   [result.Errors, result.Verrs] = ShowResults(result,'-match',matcherr);
                elseif listonly == 4 %show only errors relating to the fullvs
                   [result.Errors, result.Verrs] = ShowResults(result,'-fullv');
                else
                   [result.Errors, result.Verrs] = ShowResults(result);
                end
            else
                result.Errors = [];
            end          
            return;
        end
        if (forceread && isempty(exptlist)) || clearerrs
            fullv.ClearErrors(FullV);
        end
        d = mydir([FullV '/*FullV.mat']);
        if isempty(d)
            result.name = FullV;
            result.errs = 'No FullVs';
        else     
            e = GetExptNumber({d.name});
            [e,id] = sort(e);
            d = d(id);
            id = find(e>0);
            e = e(id);
            d = d(id);
            setid = [];
%try to avoid calling ReadExptDir if we are not doing anything            
            if ~isempty(exptlist)
                [Expts, ExptErrs] = ReadExptDir(FullV);
                [id, ia] = intersect(GetExptNumber(Expts),exptlist);
                xid = setdiff(GetExptNumber(Expts),id);
                for j = 1:length(xid)
                    EX.excluded{xid(j)} = 'Excluded on Command Line';
                end
                Expts = Expts(ia);
                EX.needcheck = id;
            elseif isfield(X,'checktimes') && ~forceread %only check times if no expts named
                eid = find(e <= length(X.checktimes));
                xid = find(X.checktimes(e(eid)) > [d(eid).datenum]); %already checked
                EX.exclude = [];
                for j = 1:length(xid)
                    EX.excluded{e(eid(xid(j)))} = sprintf('Not changed since Check on %s',datestr(X.checktimes(e(eid(xid(j))))));
                    EX.exclude(e(eid(xid(j)))) = 1;
                end
                id = find(X.checktimes(e(eid)) < [d(eid).datenum]); %not checked
                needexpts = [e(eid(id)) e(e > length(X.checktimes))];
%                d = d(eid(id)); %keep file list - report reasons for NOT
%                loading
                 if checkmains
                     if isfield(X,'exptdetails') && isfield(X.exptdetails,'checkmains')
                     else
                         needexpts = e;
                     end
                 end
                 if isempty(Expts) && ~isempty(needexpts)
                     Expts = ReadExptDir(FullV);
                     [id, ia] = intersect(GetExptNumber(Expts),needexpts);
                     Expts = Expts(ia);
                 end
                 EX.needcheck = [id];
                         
                nexpts = max([id(:)']);
                xid = setdiff(1:nexpts,[id find(EX.exclude)]);
                for j = 1:length(xid)
                    EX.excluded{xid(j)} = sprintf('No Expt for %d',xid(j));
                end
                
            elseif isempty(Expts) %either forced or no existing record of times. Do all
                [Expts, ExptErrs] = ReadExptDir(FullV);
                EX.needcheck = GetExptNumber(Expts);
            end
            if ~isempty(Expts) || forceread
                BackupFile(errfile,'print','copy');
                result = fullv.Check({d.name},Expts,EX,'nobackup', varargin{:});
                new = length(Expts);
            else %all fullv files have been checked
                new = 0;
            end
            if new == 0 && verifyerrors >0 %no new checks. Show Old Errors
                [Expts, ExptErrs] = ReadExptDir(FullV);
                if isfield(ExptErrs,'Errors') && ~isempty(ExptErrs.Errors);
                    result = AddError(result, ExptErrs);
                    fullv.SaveErrors(result);
                end
            elseif new == 0
                fprintf('Previous  Errors for %s ',FullV);
                [~,result.Verrs] = ShowResults(result);
            else
                if forceread %new list
                    args = {};
                else
                    args = {};
                end
                SaveNewErrors(FullV,result,'nobackup',ExptErrs,args{:});
            end
            if clearspk && sum(CellToMat(result,'nerr')) == 0
               d = mydir([FullV '/*spkblk*.mat']);
            end
        end
    else %char but not dir
        if ~exist(FullV) && ~exist(fileparts(FullV))
            fprintf('Cant Find Files for %s\n',FullV);
            return;
        end
        if isempty(Expts)
            Expts = ReadExptDir(fileparts(FullV));
        end
        if isempty(Expts)
            fprintf('No Expts in %s (for %s)to check\n',fileparts(FullV),FullV);
            return;
        end
        exptid = GetExptNumber(Expts);
        id = find(GetExptNumber(Expts) == GetExptNumber(FullV));
        if ~isempty(id)
            Expt = Expts{id};
        end
        fprintf('Checking %s\n',FullV);
        X = LoadFullV(FullV);
        result = fullv.Check(X,Expt,varargin{:});
        FullV = X;
        if isfield(X,'errmsg') && ~isempty(X.errmsg) %errors recorded in file
            F = X.errdata;
            if iscell(F) && length(F) == 2 && isempty(F{1}) && isstruct(F{2})
                F = F{2};
            end
            for j = 1:length(F)
                F(j).s = X.errmsg{j};
            end
            result.FullVFileErrors = F;
            if clearerrs && result.nerr == 0
                result.olderrs = F;
                FullV.olderrs = F;
                FullV.errmsg = {};
                FullV.errdata = [];
                fullv.SaveErrors(FullV,result,'clearerr'); 
                fullv.save(FullV);
            end
        end
    end
    return;
elseif isfield(FullV,'V');
        if isempty(Expts) && isempty(Expt)
            if isfield(FullV,'loadname')
                Expts = ReadExptDir(fileparts(FullV.loadname));
            end
        end
        exptid = GetExptNumber(Expts);
        id = find(GetExptNumber(Expts) == GetExptNumber(FullV));
        if ~isempty(id)
            Expt = Expts{id};
        end
elseif isfield(FullV,'errmsg')
    result.nerr = length(FullV.errmsg);      
end

result.exptno = GetExptNumber(FullV);
%DO NOT add fields to result above here - breaks references to result{j} =
%fullv.check above
%Once here, FullV should be FullV struct
result.missingtrials = [];
result.checktime = now;
result.program = 'fullv.Check';
if ~isfield(FullV,'blkend') && isfield(FullV,'blkstart')
    for j = 1:length(FullV.blkstart)
        FullV.blkend(j) = FullV.blkstart(j)+ FullV.blklen(j).*FullV.samper;
    end
    result.blkend = FullV.blkend;
end
if ~isfield(FullV,'V') || size(FullV.V,1) < 2 %may be FullVdata from cluster file or single channel from Utah
    if isfield(FullV,'blkstart')
        if ~isfield(FullV,'V')
            fprintf('No Voltages - checking blkstarts against Expt\n');
        elseif ~isfield(FullV,'t')
            fprintf('Building times list from blkstarts\n');
        end
        t = BuildFullVt(FullV);
        [a,result.missingids{1}] = FindMissingTrials(Expt, t);
        result.missingid = [];  % can't do individual probes
        if ~isempty(a)
            fprintf('missing ids%s (t = %s)\n',sprintf(' %d',result.missingids{1}),sprintf(' %.3f',result.missingids{1}));
        end
    else
        fprintf('No Voltages to check\n');
    end
    return;
end

if ~isfield(FullV,'chstd')
    FullV.chstd = std(FullV.V,[],2);
end
if checkmains
    Ex = my.load(BuildFileName(FullV,'spike2mat'),'safe');
    mains = Spike2.FindChannel(Ex,'mains');
    if isfield(mains,'times')
        result.checkmains = CheckMains(FullV, mains.times);
    end
end
result.size = size(FullV.V);
result.chstd = FullV.chstd;
%at 40Khz this is about the range spanned in 1 sec of data
result.chrange = diff(prctile(FullV.V',[0.01 99.9]));
if isfield(FullV,'loadname') %absent if being built from spkblk
    result.name = FullV.loadname;
else
    result.name = FullV.name;
end
ends = cumsum(FullV.blklen);
if isfield(FullV,'exptno')
    result.exptno = FullV.exptno;
else
    result.exptno = GetExptNumber(Expt);
end
starts = [1 1+ends(1:end-1)];
errorfile = [];
if isfield(FullV,'loadname')
    errorfile = [fileparts(FullV.loadname) '/Errors.mat'];
    BackupFile(errorfile,'copy');
end
%When thre are large artifacts, std meanV can be bigger than individual
%channels. loa
FullV = fullv.fix(FullV,'todouble');  %make sure its double
if isfield(FullV,'meanV')
    
%first find high negative correlations with meanv before adding it back in
%in case channel was blank and then had MeanV removed
    
    for j = 1:length(FullV.blklen)
        for k= 1:size(FullV.V,1)
            xc = corrcoef(FullV.V(k,starts(j):ends(j)),FullV.meanV(starts(j):ends(j)));
            meanxc(j,k) = xc(1,2);
        end
    end
    [a,b] = find(meanxc < -0.8);
    for j = 1:length(a)
        result = AddError(result,'-write',errorfile,'Expt %d. Some Channels (%s) close to MeanV. Blocsk %s)',...
            FullV.exptno,sprintf(' %d',unique(a)),sprintf(' %d',unique(b)));
    end
       
    
    crit = std(FullV.meanV)./10;
%add back in hte mean voltage for channels where variance is bigger
%in case the channel is actually blank
    if isfield(FullV,'meangain')
    tsd = std(FullV.V,[],2);
    for j = 1:size(FullV.V,1)
        if tsd(j) > crit
            FullV.V(j,:) = FullV.V(j,:) + FullV.meanV * FullV.meangain(j);
        end
    end
    end

else 
    crit = 0;
end

if isfield(FullV,'loadtype') && strcmp(FullV.loadtype,'double')
    cprintf('red', 'Expt %d Already double on disk\n',FullV.exptno);
    vmax = max(abs(FullV.V(:)));
    if fix.int
        fix.save = 1;
    end
    if ~isfield(FullV,'intscale')
        FullV.intscale(1) = vmax;
        FullV.intscale(2) = double(intmax('int16')-5);
    end
    if  vmax < FullV.intscale(1)/2 %rescaled
        cprintf('red', 'Expt %d Badly scaled Max %f, should be %.4f\n',FullV.exptno,vmax,FullV.intscale(1));
        if fix.scale || fix.int
            rescale = FullV.intscale(1)./vmax;
            FullV.V = FullV.V .* rescale;
            FullV = fullv.ClearErrors(FullV,{'flat' 'missingprobe'});
            fix.save = 1;
        end
    end
end


meanxc = [];
for j = 1:length(FullV.blklen)
    sds(j,:) = std(FullV.V(:,starts(j):ends(j))');
end

%first find high correlations with meanv
[a,b] = find(abs(meanxc) > 0.8);
for j = 1:length(a)
    result = AddError(result,'-write',errorfile,'Expt %d. Some Channels (%s) close to MeanV. Blocsk %s)',...
    FullV.exptno,sprintf(' %d',unique(a)),sprintf(' %d',unique(b)));    
end
   
%prctile sds 80 should ensuer that this is the value of a read channle.
%Channels with amplitude < 1/10 of this are suspicious
if crit == 0
    crit = prctile(sds(:),80)./10;
end

std(sds);
err(1) = 0;
if isfield(Expt,'Trials')
    [result.missingtrials, result.missingid, goodid] = FindMissingTrials(Expt,FullV.t.*10000);
    for j = 1:length(Expt.Trials)
        trials(j,1) = Expt.Trials(j).Start(1)./10000;
        trials(j,2) = Expt.Trials(j).End(end)./10000;
        trialids(j) = Expt.Trials(j).id;
    end
    if ~isempty(result.missingid)
        result = AddError(result,'-write',errorfile,'Expt%d missing Trial ids %s (%s)',FullV.exptno,...
            sprintf(' %d',result.missingid),sprintf(' %.1f',trials(result.missingtrials,1)));
    end
    for j = 1:size(goodid,1)
    end
else
    trials = [];
end

[a,b] = find(sum(sds(:,1:16) < crit,2) > 14);
if ~isempty(a)
    result = AddError(result,'-write',errorfile,'Missing probes 1:16 in Expt%d Blocks %s',FullV.exptno,sprintf(' %d',a));    
    err(1)  = 1;
    if ~isempty(Expt)
        result.missinggroup{1} = FindTrialsInBlocks(FullV,a, trials);
        result.missingids{1} = trialids(result.missinggroup{1});
    end
    fix.save = 0;
end
if size(sds,2) < 17
    a = 1:size(sds,1);
else
    [a,b] = find(sum(sds(:,17:24) < crit,2) > 7);
end
if ~isempty(a)
    result = AddError(result,'-write',errorfile,'Missing probes 17:24 in Expt%d Blocks %s',FullV.exptno,sprintf(' %d',a));
    result.missinggroup{2} = FindTrialsInBlocks(FullV,a, trials);
    result.missingids{2} = trialids(result.missinggroup{2});
    err(2) = 1;
    fix.save = 0;
end

if sum(err) == 0 && sum(sds(:) < crit)
    [a,b] = find(sds < crit); %prob missing
    bad = unique(b); %blocks with errors
    for j = 1:length(bad)
        id = find(b == bad(j));
        result = AddError(result,'-write',errorfile,'FullV Flat in Blocks %s in Expt%d Probe %d',sprintf(' %d',a(id)),FullV.exptno,bad(j));
    end
    [a,b] = find(sds < crit); %low amp
    weak = unique(b); %blocks with errors
    weak = setdiff(weak,bad);
    for j = 1:length(weak)
        id = find(b == weak(j));
        result = AddError(result,'-write',errorfile,'FullV small in Blocks %s in Expt%d Probe %d',sprintf(' %d',a(id)),FullV.exptno,bad(j));
    end
end

if isfield(result,'errs')
    result.nerr = length(result.errs);
else
    result.nerr = 0;
end
if fix.save
    fullv.save(FullV,'verbose');
end
fullv.SaveErrors(FullV,result,sargs{:}); %also record any in thez file

function SaveNewErrors(name, res, varargin)

E.errmsg = {};
E.errdata = [];
E.loadname = [name '/Errors.mat'];


%fields to copy from res
cpf = {'chstd' 'size' 'checkmains' 'chrange' 'checktime' 'nerr'};
for j = 1:length(res)
    if isfield(res{j},'errs')
        E = AddError(E,res{j});
    end
    if isfield(res{j},'exptno') %res{j} may be empty
    e = res{j}.exptno;
    if e > 0
        E.exptdetails(e).exptno = e;
        if isfield(res{j},'checktime')
            E.checktimes(e) = res{j}.checktime;
        end
        for k = 1:length(cpf)
            f = cpf{k};
            if isfield(res{j},f)
                E.exptdetails(e).(f) = res{j}.(f);
            end
        end
    end
    end
end
for j = 1:length(varargin)
    if isfield(varargin{j},'Errors')
        E = AddError(E, varargin{j});
    end
end
X = fullv.SaveErrors(E,varargin{:});


function tlist = FindTrialsInBlocks(FullV, bid, trials)

tlist = [];
if isempty(trials)
    return;
end
for j = 1:length(bid)
    tid = find(trials(:,1) > FullV.blkstart(bid(j)) & trials(:,2) < FullV.blkend(bid(j)));
    tlist = [tlist tid'];
end




function [errs, verrs] = ShowResults(res, varargin)
strargs = cell2cellstr(varargin);
showmains = 1;

%find errors that suggest there are problems with FUllV Files that
%might be fixed with new spkblks
errtype = [];
errs ={};
if iscell(res)
    return;
end
    
res.Errors = FixErrors(res.Errors);
errcount = zeros(1,3);
for j = 1:length(res.Errors)    
    s = res.Errors(j).s;
    if isfield(res.Errors,'progname') 
        if sum(strcmp(res.Errors(j).progname,{'APlaySpkFile' 'expt.Check' 'SetExptRC'}))
            errtype(j) = 0;
        elseif sum(strcmp(res.Errors(j).progname,{'fullv.Check' 'LoadFullV' 'BuildFullV'}))
            errtype(j) = 1;
            if strfind(s,'Excluded Time Range')
                errtype(j) = 4;
                id = regexp(s,'Chans at [0-9]+');
                if ~isempty(id)
                    t = sscanf(s(id(1)+9:end),'%f-%f');
                    if diff(t) < 4
                        errtype(j) = -1;
                    end
                end
            elseif  regexp(s,'Missing [0-9]+ Trials for')
                errtpe = 5;
                t = sscanf(s(9:end),'%f');                
                if t < 2 %just one missing trial - not a spkblk error
                    errtype(j) = -5;
                end
            elseif strfind(s,'missing Trial ids')
                id = strfind(s,'missing Trial ids');
                t = sscanf(s(id(1)+18:end),'%f');
                if length(t) < 2 %just one missing trial - not a spkblk error
                    errtype(j) = -1;
                end
            elseif regexp(s,'Trial [0-9]+.* missing probes')
                id = regexp(s,'Trial [0-9]+.* missing probes');
                errcount(1) = errcount(1) +1;
                if errcount(1) > 2
                    errtype(j) = 6;
                else
                    errtype(j) = -6;
                end
            end
        else
            errtype(j) = 2;
        end
    else
        errtype(j) = 3;
    end
end

verrs = find(errtype > 0);
if sum(strcmp('-fullv',strargs))
    if ~isempty(verrs)
        errs = ShowErrors(res, '-match', {'FullV Flat' 'Cannot Read' 'Excluded Time Range' 'blocks end\(n\)' 'missing Trial' 'Missing probes'});
        [a,b,c] = MakeProbeIndex(fileparts(res.name));
        for j = 1:length(b)
            s = CheckExceptions(b{j}.exc,'-silent');
            errs = AddError(errs,'-silent','-prog','Spike2',s);
        end
    else
        errs = {};
    end
    showmains = 0;
else
    errs = ShowErrors(res, varargin{:});
end
if isfield(res,'exptdetails') && sum(strcmp('-match',strargs)) == 0 && showmains
    if isfield(res.exptdetails,'checkmains')
        x = cat(1,res.exptdetails.checkmains);
        a = cat(1,x.amps);
        c = mean(a);
        amax = max(a);
        if isfield(x,'meanamp')
            m = cat(1,x.meanamp);
            b = max(m);
            fprintf('Mains Noise Max(mean) Amp: %.3f(%.3f),%.3f\n',amax,c,b);
        end
    end
    if isfield(res.exptdetails,'chstd')
        ne = 0;
        for j = 1:length(res.exptdetails)
            if ~isempty(res.exptdetails(j).chstd)
                ne = ne+1;
                nch(ne) = length(res.exptdetails(j).chstd);
                V(ne,1:nch(ne)) = res.exptdetails(j).chstd;
            end
        end
        a = mean(V,1);
        fprintf('Signal Amplitude %.2f  +- %.2f: %.2f - %.2f\n',mean(a),std(V(:)),min(a),max(a));
        if length(unique(nch)) > 1
            cprintf('red','!!!!!!!!! Varaible Channel Count !!!!!!!!!!!!!!!\n')
            cprintf('red','%d ',nch);
            fprintf('\n');
        end
    end
end


function X = CheckMains(FullV, tmains)

tic;
period = round(0.01662./FullV.samper);
tv = 1:(period+10);
vt = FullV.t;
ts(length(vt)) = 0;
d = tmains - FullV.blkstart(1);
id = find(d > 0);
        t = tmains(id(1));
        tid = find(vt >= t,1);
        ti = tid(1); % first sampale
        x = ti:ti+tv;
        ts(ti+tv-1) = tv;
        ts(1:ti-1) = tv(end-ti+2:end);
        mv = FullV.V(:,ti+tv-1);
        if isfield(FullV,'meanV')
            mmv = FullV.meanV(ti+tv-1);
        end

        mt(1) = ti;
        j = 2;
        while t < vt(end) && ti+2*period+100 < length(vt) && j < length(id)
            tid = find(vt(ti+period-5:ti+period+5) >= tmains(id(j)),1);
            if isempty(tid)
                tid = find(vt(ti+period-5:end) >= tmains(id(j)),1);
            end
            ti = tid(1)+ti+period-5;
            x = ti:ti+tv;
            ts(ti+tv-1) = tv;
            mv = mv+FullV.V(:,ti+tv-1);
            if isfield(FullV,'meanV')
                mmv = mmv + FullV.meanV(ti+tv-1);
            end
            mt(j) = ti;
            j = j+1;
        end
        ts(ti+tv(end):end) = 1;
        mv = mv./j;
        mv = mv-(repmat(mean(mv,2),1,size(mv,2)));
X.dur = toc;
if isfield(FullV,'chstd')
    X.amps = std(mv,[],2)./FullV.chstd;
else
    X.amps = std(mv,[],2)./std(FullV.V,[],2);
end
if isfield(FullV,'meanV')
    mmv = mmv./j;
    X.meanamp = std(mmv)./std(FullV.meanV);
end
ff = abs(fft(mv(1,:)));
