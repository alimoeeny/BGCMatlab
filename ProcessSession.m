function varargout = ProcessSession(dirname, varargin)
%ProcessSession(dir) calls all the other programs needed to build Cluster
%files etc. Apart from Parallelization, It should make most decisions automatically,
%
%ProcessSession(dirname,'parallel') will use parfor loops for the steps that
%benefit from this - making spkblkfile, automatic clustering
%
%If run on a PC, will call Spike2 to build spkblk files
%On Mac Needs Spike2 makemat to have run first.
% Then calls BuildAllFullV, RunAllVPcs as needed
%
%
%ProcessSession(dirname,'refcut') runs a refcluster cut on all the Fullvs
%          Does this automatically if there is a RefClusters file
%          ...,'norefcut') prevents refcuts from being made.
%
%                ...,'noauto') does not do automatic cluster
%                ...,'skipfullv') do not make fullv files (assumes they exist)
%                ...,'scan') passes this arg to BuildAllFullV, to run when Spike2 is still
%              making fullspk files. Not needed on PCs
args = {};
strargs = cell2cellstr(varargin);

cleanspkblk = 1; %now the default. Remove spkblk files once FullV are done and no errors
if sum(strcmp('cleanup',strargs))
    cleanspkblk = 1;
end

if iscellstr(dirname)
    for j = 1:length(dirname)
        fprintf('\n\n###################################\nProcessSession For %s:\n',dirname{j});
        log{j} = ProcessSession(dirname{j},varargin{:});
    end
    varargout{1} = log;
    return;
elseif iscell(dirname)
    for j = 1:length(dirname)
        ProcessSession(dirname{j},varargin{:});
    end
    return;
end

if isfield(dirname, 'expstate') %A log
    SummarizeLog(dirname,varargin{:});
    return;
elseif isdir(dirname)
    sessionlog = [dirname '/SessionLog.mat'];
    if exist(sessionlog,'file')
        Session= my.load(sessionlog,'safe');
        Session.Logs{end+1}.user = GetUserName;
    else
        Session.Logs{1}.user = GetUserName;
        Session.myload.name = sessionlog;
    end
    if sum(strcmp('showlog',strargs))
        for j = 1:length(Session.Logs)-1
            SummarizeLog(Session.Logs{j});
        end
        varargout{1} = Session.Logs(1:end-1);
        return;
    end
    Session.Logs{end}.name = dirname;
    Session.Logs{end}.args = strargs;
    Session.Logs{end}.events.start = now;
    d = dir([dirname '/RefClusters.mat']);
    
    if ~isempty(d) && sum(strncmp('noref',strargs,5)) == 0
        fprintf('Found %s. Will build Refclusters if needed\n', [dirname '/RefClusters.mat']);
        fprintf('(Disable this with ''noref'' argument)\n');
        strargs = {strargs{:} 'refcut'};
    end
    fprintf('Checking for existing files in %s\n',dirname);
    [need, matfiles] = CheckDirInit(dirname,strargs{:});
    CheckProbeIndex(dirname, matfiles);
    Session.Logs{end}.events.dirdone = now;
    
    Session.Logs{end}.need = need;
    if sum(strcmp('parallel',strargs))
        args = {args{:} 'parallel'};
    end
%Build all matfiles first, and run ReadExptDir, since these
%steps often generate errors. Rest shoudl be good to run overnight
    if need.matfile || (need.spkblk && need.fullv)
        if ispc
            fprintf('Building SpkBlk files with Spike2\n');
            spk2res = Spike2.MakeSpkBlk(dirname, 'mklfp',args{:});
        else
            fprintf('Need to Run on a PC, or run MakeMat in SPike2 first\n');
            return;
        end
    end
    Session.Logs{end}.events.spkblkdone = now;
    sessionstart = now;
%Now .mat files should exist for all Expts, Convert these wiht
%APlaySpkfile, and build a list of Expts
    [Expts, details] = ReadExptDir(dirname);
    Session.Logs{end}.expstate = details;
    Session.Logs{end}.events.readexptsdone = now;
    exptnos = GetExptNumber(Expts);
    if sum(strcmp('scan',strargs))
        bargs = {'scan'};
    else
        bargs = {};
    end
    
%Build FullV Files for al    
    if sum(strcmp('skipfullv',strargs))
        fprintf('Skipping FullV Build\n');
    elseif sum(strcmp('checkonly',strargs))
        if sum(details.needfullv) ==0
            fprintf('FullV are up to date\n');
        else
            fprintf('Need FullV For Expts:%s\n',sprintf(' %d',details.needfullv));
        end
    elseif sum(details.needfullv) %ReadExptDir says these are missing
        fprintf('Indexing contexts of spkblk files');
        [probes, errs, probeerrs] = MakeProbeIndex(dirname);
        idir = dir([dirname '/' 'probes.mat']);
        e = [];
        id = [];
        if isfield(matfiles,'spkblk')
            for j = 1:length(matfiles.spkblk)
                e(j) = GetExptNumber(matfiles.spkblk{j});
            end
            id = find(matfiles.spkblktime > idir.datenum);
        end
%this will detect any files where no spkblkfiles were made last time, but not if they were updated
        missing = setdiff(unique(e),unique([probes.suffix]));
        if sum(missing > 0) || ~isempty(id) && 0
            fprintf('Rebuilding index of SpkBlk files for Expts %s\n',sprintf('%d ',missing));
            MakeProbeIndex(dirname,'reindex');
            Session.Logs{end}.events.probeindexdone = now;
        end            
        if sum(missing >0) && isempty(matfile.spkblk)
            if ispc
                fprintf('No spkblk files. Attempting to rebuild\n');
                eid = find(details.needfullv);
                spk2res = Spike2.MakeSpkBlk(dirname, 'mklfp','expts',eid);
                MakeProbeIndex(dirname,'reindex');
            else
                fprintf('No spkblk files. Need to run on a PC\n');
            end
           
        end
        fprintf('Building FullV Files\n');
        ts = now;
        BuildAllFullV(dirname,'nocut',bargs{:});
        Session.Logs{end}.events.fullvdone = now;
        Session.Logs{end}.dur.fullv = mytoc(ts);
        
    else
        fprintf('All FullV Files (Expts%d - %d) Exist\n',min(exptnos),max(exptnos));
    end

    if ~sum(strcmp('doauto',strargs)) %doauto forces check of autoclusters
        details.needautoclusters(details.needclusters==0) = 0;
    end
    if sum(strcmp('noauto',strargs))
        details.needautoclusters = 0;
    end
    if sum(strcmp('checkonly',strargs))
        fprintf('Need Auto Clusters for Expts:%s\n',sprintf(' %d',details.needautoclusters));
    elseif sum(details.needautoclusters)
        e = intersect(exptnos,find(details.needautoclusters));
        ts = now;
        RunAllVPcs(dirname,'expts',e,'autocutall','savespikes',args{:}); 
        Session.Logs{end}.events.autocutdone = now;
        Session.Logs{end}.dur.autocut = mytoc(ts);
    elseif sum(strcmp('checkclusters',strargs))
        fprintf('Checking ClusterTimes Files\n');
        A = GetArrayConfig(dirname,'guess');
        nprobes = length(A.id);
        for j = exptnos(:)'
            Clusters = {};
            cname = BuildFileName(dirname,'autocluster',j);
            load(cname);
            if isempty(Clusters) || length(Clusters) < nprobes
                fprintf('Only %d probes in %s\n',length(Clusters),cname);
                RunAllVPcs(dirname,'expts',j,'autocutall','savespikes',args{:});
            end
        end
    else
        fprintf('AutoClusters Finished\n');
    end
    if sum(strcmp('checkonly',strargs))
        fprintf('Need Manual/Ref Clusters for Expts:%s\n',sprintf(' %d',details.needclusters));
    end
    if sum(strcmp('refcut',strargs))
        e = intersect(exptnos,find(details.needclusters));
        if isempty(e)
            fprintf('ClusterTimes Files all present. To force remake of Refcuts use:\n');
            fprintf('RunAllVPcs(%s,''expts'',explist,''refclusterforce'',''savespikes'',''parallel'')\n',dirname);            
        else
            Array = GetArrayConfig(dirname,'guess');
            fprintf('Running RefCut');
            ts = now;
            RunAllVPcs(dirname,'expts',e,'refcut','savespikes',args{:});
            Session.Logs{end}.events.refcutdone = now;
            Session.Logs{end}.dur.refcut = mytoc(ts);
            Session.Logs{end}.expstate.needrefcut = e;
        end
    end
    if sum(strcmp('eckerfit',strargs))
        autodir = [dirname '/Ecker'];
        if ~exist(autodir)
            mkdir(autodir);
            mkdir([autodir '/Spikes']);
            eid = exptnos;
        else
            d = dir([autodir '/*AutoClusterTimes.mat']);
            if isempty(d)
                eid = exptnos;
            else
                eid = setdiff(exptnos, GetExptNumber(d));
                lastdate =  max([d.datenum]);
            end
        end            
        ts = now;
        if ~isempty(eid)
            RunAllVPcs(dirname,'expts',eid,'autocutmode','ecker','autocutall','savespikes',args{:});
        end
        celllist = [autodir '/AutoCellList.mat'];
        d = dir(celllist);
        Session.Logs{end}.events.eckerdone = now;
        Session.Logs{end}.dur.ecker = mytoc(ts);
        ts = now;
        if isempty(d) || d.datenum < lastdate %need to build/rebuild celllist
            fprintf('Building Automatic CellList\n');
            x = PlotClusters(dirname,'loadecker','autolistsave');
            PC.PlotClusters(x.toplevel,[],'close');
            close(x.toplevel);
        else
            fprintf('CellList %s is up to date\n',celllist);
        end
        Session.Logs{end}.events.autolistdone = now;
        Session.Logs{end}.dur.autolist = mytoc(ts);
    end
    
%If no refcut, then Copy AutoClusterTimes to ClusterTimes, if ClusterTimes does not exist    
    if sum(strcmp('refcut',strargs)) == 0
        InitClusters(dirname);
    end
    ts = now;
    Session.Logs{end}.events.startcheck = ts;
    CD  = [];
    if sum(details.needfullv) ==0 && sum(details.needautoclusters) ==0 && sum(details.needclusters) ==0
        if sum(strcmp('nocheck',strargs)) == 0
            fprintf('Checking Integrity of ClusterTimes Files in %s\n',dirname);
            CD = CheckExptClusters(dirname,Expts, 'nocells','quiet');
        end
    else
        CD = CheckExptClusters(dirname,Expts, 'nocells','quiet');
    end
    if ~isempty(CD)
        id = find(CD.nspks ==0);
        if isempty(id)
            fprintf('All Expts/Probes have Clusters Defined\n');
        else
            fprintf('Exps with missing Spikes:\n');
            [e,p] = ind2sub(size(CD.nspks),id);
            for j = 1:length(e)
                fprintf('Expt%dP%d %d/%d Events\n',e(j),p(j),CD.nspks(e(j),p(j)),CD.nevents(e(j),p(j)));
            end
        end
    end
    try
        ListExpts(dirname,'silent');
    end
    
    
    a = BuildRFFits(dirname);
    if isfield(a,'savedate') && a.savedate > ts %new RF data
        [mdir, sname] = fileparts(dirname);
        if regexp(sname,'M[0-9]+')
            BuildRFFits([mdir '/M*'],'save');
        end
    end
    
    
    ShowErrors(details,'severe');
    ShowErrors(CD,'severe');
    if isfield(CD,'program')
        Session.Logs{end}.ClusterCheck = CD;
    end
    X = fullv.Check(dirname);
    go = 1;
    s = 'unknown case';
    spkcheck.spkerr = 0;
    if ~isempty(X.Verrs)
        go = 0;
        s = sprintf('Check FullV had %d Errors',length(X.Verrs));
        spkcheck.spkerr = 3;
    end
    if need.fullv %will have build new fullv files. Build list again
         [need, matfiles] = CheckDirInit(dirname,strargs{:});
    end
    
    if  ~isempty(matfiles.spkblk);
        spkcheck.spkerr = 1;
        try
        if ~isfield(X,'exptdetails') || ~isfield(X.exptdetails,'chstd')
            go = 0;
            s = sprintf('FullV Check did not include Voltage SD');
        else
           spke = unique(GetExptNumber(matfiles.spkblk));
           fulle = GetExptNumber(matfiles.fullv);
           missing = setdiff(spke,fulle);
           if ~isempty(missing)
               s = sprintf('Missing FullV for Expts %s',sprintf('%d ',missing));  
               go = 0;
           end
        end
        if go
           chsd = cat(2,X.exptdetails.chstd);
           crit = mean(chsd(:)) - std(chsd(:)) .* 3;
           if crit < 0.01
               crit = 0.01; %std < 0.01V is suspicious
           end
           [a,b] = find(chsd<crit);
           if isempty(a)
               s = sprintf('FullV files All Correct (%d)',length(fulle));
               spkcheck.spkerr = 0;
               if cleanspkblk
                   for j = 1:length(matfiles.spkblk)
                       fprintf('Deleting %s\n',matfiles.spkblk{j});
                       delete([dirname '/' matfiles.spkblk{j}]);
                   end
                   s = sprintf('Deleted All (%d) spkblk files',length(matafiles.spkblk));
               end
           else               
               spkcheck.spkerr = 2;
               s = sprintf('FullV Channels with Low amplitude:',length(a));
               fprintf('%s (mean is %.4f)\n',s,mean(chsd));
               for j = 1:length(a)
                   fprintf('E%dP%d %.4f\n',a(j),b(j),chsd(a(j),b(j)));
               end
               spkcheck.lowamp = [a b];
           end
        end
        fprintf('%d spkblk Files: %s\n',length(matfiles.spkblk),s);
        spkcheck.result = s;
        Session.Logs{end}.spkblkcheck = spkcheck; 
        catch ME
            CheckExceptions(ME);
        end
    else
        spkcheck.result = 'No Spkblk files';
    end
    Session.Logs{end}.spkblkcheck = spkcheck;
    
    
    Session.Logs{end}.dur.check = mytoc(ts);
    Session.Logs{end}.events.finish = now;
    my.save(Session,'-safe','-verbose');
    varargout{1} = Session.Logs{end};
end

function InitClusters(name)
%Copy AutoClusters over to Clusters so that PlotClusters can be run without
%runauto
if ~isdir(name)
    fprintf('%s is not a folder\n',name);
end

d = mydir([name '/*AutoClusterTimes.mat']);
for j = 1:length(d)
    ename = strrep(d(j).name,'AutoCluster','Cluster');
    if ~exist(ename)
        fprintf('%s copied to %s\n',d(j).name,ename);
        copyfile(d(j).name,ename);
        xname = strrep(d(j).name,'AutoClusterTimes','AutoClusterTimesDetails');
        ename = strrep(xname,'AutoCluster','Cluster');
        copyfile(xname,ename);
        fprintf('%s copied to %s\n',xname,ename);
    end
end


function [need, files] = CheckDirInit(name, varargin)
%check that initial files have been built if necessary

if sum(strcmp('checkonly',varargin))
    fullcheck = 1;
else
    fullcheck = 0;
end
need.matfile = 0;
need.lfp = 0;
need.spkblk = 0;
need.fullv = 0;
files = [];

d = dir([name '/*.mat']);
if isempty(d)
    need.lfp = 1;
    need.spkblk = 1;
    need.matfile = 1;
    need.fullv = 1;
    return;
end

id = find(CellToMat(strfind({d.name},'FullV.mat')));
if isempty(id)
    need.fullv = 1;
    files.fullv = {};
else
    files.fullv = {d(id).name};
end
e = [];
id = find(CellToMat(regexp({d.name},'.*spkblk.*.mat')));
if isempty(id) && need.fullv
    need.spkblk = 1;
else
    files.spkblk = {d(id).name};
    files.spkblktime = [d(id).datenum];
    for j = 1:length(id)
        e(j) = GetExptNumber(d(id(j)).name);
    end
    files.gotspkblk = unique(e);
end
if need.spkblk ==0 || fullcheck
    sm = dir([name '/*.smrx']);
    e = [];
    for j = 1:length(sm)
        e(j) = GetExptNumber(sm(j).name);
    end
    files.smrexpts = unique(e);
    missing = setdiff(files.smrexpts,files.gotspkblk);
    if ~isempty(missing)
        files.needspkblk = missing;
        need.spkblk = 1;
        if fullcheck
            fprintf('Need SpkBlk for:');
            for j = 1:length(files.needspkblk)
                fprintf(' %d',files.needspkblk(j));
            end
            fprintf('\n');
        end
    end
end

id = find(CellToMat(regexp({d.name},'.lfp.mat$')));
if isempty(id)
    need.lfp = 1;
else
    files.lfp = {d(id).name};
end

function CheckProbeIndex(dirname, X)
%Check to see if any spkblk files are not in probe index

pfile = [dirname '/probes.mat'];
if ~exist(pfile,'file')
    return;
end
if isempty(X.spkblk)
    return;
end
P = load(pfile);
np = max([P.probes.probe]);
for j = 1:length(X.spkblk)
    e(j) = GetExptNumber(X.spkblk{j});
    if regexp(X.spkblk{j},'[0-9]+A.[0-9]+spkblk')
        suff(j) = 2;
    else
        suff(j) = 1;
    end
end
expts = unique(e);
for j = 1:length(expts)
    id = find([P.probes.suffix] == expts(j));
    p = unique([P.probes(id).probe]);
    if max(p) < np
        eid = find(e == expts(j));
        suffs = unique(suff(eid));
        if length(suffs) ==1
            fprintf('Expt%d Probe index missing some spkblk files. Only have Suff %d\n',expts(j),suffs);
        end
    end
end

function SummarizeLog(Log, varargin)

strargs = cell2cellstr(varargin);
if isfield(Log,'name')
    fprintf('%s: ',Log.name);
end
if sum(strcmp('spkblk',strargs))
    if isfield(Log,'spkblkcheck')
        fprintf('%d FullV Errors. ',Log.spkblkcheck.spkerr);
        fprintf('%s',Log.spkblkcheck.result);
    end
    fprintf('\n');
    return;
end

E = Log.events;
D = Log.dur;
fprintf('%s to %s: %.3f hours\n',datestr(E.start),datestr(E.finish),(E.finish-E.start).*24);
if isfield(D,'fullv');
    fprintf('Made %d Fullv Files. %.3f hours\n',sum(Log.expstate.needfullv),D.fullv./(60*60));
end
if isfield(D,'autocut');
    fprintf('Made %d AutoCluster Files. %.3f hours\n',sum(Log.expstate.needautoclusters),D.autocut./(60*60));
end
if isfield(D,'refcut')
    fprintf('Made %d Refcluster Files. %.3f hours\n',sum(Log.expstate.needautoclusters),D.refcluster./(60*60));
end

f = fields(E);
for j = 1:length(f)
    dates(j) = E.(f{j});
end
[~,id] = sort(dates);
for j = 1:length(id)
    fprintf('%s %.2f hrs %s\n',datestr(dates(id(j))),mytoc(dates(id(1)),dates(id(j)))./3600,f{id(j)});
end
