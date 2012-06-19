function res = RunAllVPcs(name, varargin)
%RunAllVPcs(dir, .....
%Reads all FullV Matlab files in a directory and calls AllVPcs

expts = [];
bysuffix = 0;
scanfiles = 0;
recalc = 0;
sizecheck = 0;  %set to max safe bytes for FullV to control if FullV is kept in memory
runautocut = 1;
%memsize = CheckPhysicalMemory;
%sizecheck = (memsize * 1e6) .* 0.8;
X.spkrate = 50;
template = [];
forcerun = 0;
parallel = 0;
args = {};
addargs = {};
checktype = 'reclassify';

if ischar(name)
    logfile = [name '/' 'RunAll.log'];
    logfid = fopen(logfile,'a');
elseif isfield(name,'prefix')
    logfile = [name.prefix '/' 'RunAll.log'];
    logfid = -1;
else
    logfid = -1;
end
X.logfid = logfid;


j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
    elseif strncmpi(varargin{j},'args',4)
        j = j+1;
        addargs = varargin{j};
    elseif strncmpi(varargin{j},'checktype',9)
        j = j+1;
        checktype = varargin{j};
    elseif strncmpi(varargin{j},'expts',3)
        j = j+1;
        expts = varargin{j};
    elseif strncmpi(varargin{j},'nocut',4)
        runautocut = 0;
        args = {args{:} 'nocut'};
    elseif strncmpi(varargin{j},'paralell',4)
        parallel = 1;
    elseif strncmpi(varargin{j},'run',3)
        forcerun = 1;
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end
res.starts(1) = now;
warning('off','stats:gmdistribution:FailedToConverge');
warning('off','stats:gmdistribution:MaxIterations');

if iscell(name) &&  length(name) > 1
     name = CellToStruct(name);
end

     
if isstruct(name)
    if isfield(name,'name') && isfield(name,'probes') && isfield(name,'args');
        if isfield(name,'overlapn') && forcerun == 0
            res = CheckReclassify(name,'reclassify');
        else
            res = RunCommandList(name, varargin{:});
        end
    elseif isfield(name,'overlapn') && forcerun == 0
            res = CheckReclassify(name,'reclassify');
    elseif isfield(name,'acts')
        if ~isfield(name.acts,'name')
            for j = 1:length(name.acts)
                name.acts(j).name = sprintf('%s/%s',name.prefix,name.name{name.acts(j).exptid});
            end
        else
            for j = 1:length(name.acts)
                if strncmp(name.acts(j).name,'Expt',4)  %no directory prefix
                    name.acts(j).name = sprintf('%s/%s',name.prefix,name.acts(j).name);
                end
            end
        end
            res = RunCommandList(name.acts,'log',logfile,varargin{:});
    elseif isfield(name,'exptlist')
        if isfield(name,'args') &&length(addargs)
            name.args = {name.args{:} addargs{:}};
        end
        name = CheckExptList(name);
        res = CheckReclassify(name, checktype);
    end
    return;
 elseif isdir(name)
    path = name;
else
[path, fname] = fileparts(name);
end

d = dir(name);
ns = 1;
res.prefix = name;
for j = 1:length(d)
        sid = regexp(d(j).name,'Expt[0-9,a]*FullV');
        X.logfid = logfid;
        if length(sid)
            ex = sscanf(d(j).name(sid(1)+4:end),'%d');
            if isempty(expts) || ismember(ex,expts)
                fprintf(logfid,'Processing %s %s\n',d(j).name,datestr(now));
                name = [path '/' d(j).name];
                if parallel == 0
                    res.cls{ns} =  ProcessSuffix(name, ex, args{:});
                else
                    res.cls{ns} = [];
                end
 %              fclose all;
               res.exptlist(ns) = GetExptNumber(d(j).name);
               res.name{ns} = d(j).name;
               filenames{ns} = name;
               ns = ns+1;
            end
        end
        
end

if parallel
    parfor  (j = 1:ns-1)
        cls{j} =  ProcessSuffix(filenames{j}, res.exptlist(j), args{:});
    end
    res.cls = cls;
end

 res.end = now;
 res.args = args;
 if logfid > 0
     fclose(logfid);
 end
 
 
function res = ProcessSuffix(path, ex, varargin)
template = [];
nocut = 0;
res = {};
j = 1;
checkexpts = 0;
makesmall = 1;

args = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'check',5)
        checkexpts = 1;
        if strncmpi(varargin{j},'checkandbuild',8)
            checkexpts = 2;
        end
    elseif strncmpi(varargin{j},'nocut',5)
        nocut = 1;
    elseif strncmpi(varargin{j},'template',6)
        j = j+1;
        template = varargin{j};
        cargn = j;
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end

res.fiename = path;
outname = path;
if exist(outname,'file')
    ts = now;
    fprintf('%s Exits\n',outname);
    d = dir(outname);
else
    fprintf('Cant Read %s\n',outname);
    return;
end
ts = now;
res = AllVPcs(outname,args{:});
if ~iscell(res)
    res.runtime = mytoc(ts);
else

    toplevel = 0;
    for j = 1:length(res)
        if isfield(res{j},'toplevel')
            toplevel = res{j}.toplevel;
        end
        if makesmall
            res{j} = rmfields(res{j},'t','Evec','pcs');
            if isfield(res{j},'cluster')
            res{j}.cluster = rmfields(res{j}.cluster,'times');
            end
        end
    end
    if toplevel > 0 && isappdata(toplevel,'Vall')
        rmappdata(toplevel,'Vall');
    end
end



function Ex = LoadExpt(name, ei, rebuild, logfid)
        Ex = [];
        smrname = regexprep(name,'lem/M([0-9]*)','$0/lemM$1');
        exfile = [smrname '.' num2str(ei) 'idx.mat'];
        matfile = [smrname '.' num2str(ei) '.mat'];
        if exist(matfile) && (~exist(exfile,'file') || rebuild)
            PrintMsg(logfid,'Building %s\n',exfile);
            APlaySpkFile(matfile,'bysuffix','noerrs');
        end
        if exist(exfile,'file')
            fprintf('Loading %s\n',exfile);
            load(exfile);
            id = find(Expt.Trials.Result == 1);
            for t = length(id):-1:1
                Ex.Trials(t).Start = Expt.Trials.Start(id(t));
                Ex.Trials(t).End = Expt.Trials.End(id(t));
                Ex.Trials(t).ed = Expt.Trials.ed(id(t));
                Ex.Trials(t).id = Expt.Trials.id(id(t));
                Ex.Trials(t).Trial = id(t);
            end
            if isempty(ExptList)
                Ex.Header.expname = 'None';
            else
                Ex.Header.expname = ExptList(end).expname;
            end
            Ex.Header.name = exfile;
            Ex.Stimvals.ed = mean(Expt.Trials.ed(id));
        end
        
    function res = RunCommandList(X, varargin)
     
        checkmode = 0;
        logfile = [];
        addargs = {};
        j = 1;
        probes = unique([X.probes]);
        expts = unique([X.exptid]);
        while j <= length(varargin)
            if strncmpi(varargin{j},'check',5);
                checkmode = 1;
            elseif strncmpi(varargin{j},'args',4);
                j = j+1;
                addargs = varargin{j};
            elseif strncmpi(varargin{j},'expts',5);
                j = j+1;
                expts = varargin{j};
            elseif strncmpi(varargin{j},'probes',5);
                j = j+1;
                probes = varargin{j};
            elseif strncmpi(varargin{j},'log',3);
                j = j+1;
                logfile = varargin{j};
            end
            j = j+1;
        end

        rmid = zeros(size(X));
        for j = 1:length(X)
            allnames{j} = X(j).name;
            if isempty(X(j).args)
                rmid(j) = 1;
            end
        end
        id = setdiff(1:length(X),find(rmid));
        X = X(id);
        allnames = allnames(id);
        id = find(ismember([X.exptid],expts) & ismember([X.probes],probes));
        X = X(id);
        allnames = allnames(id);
        
    names = unique(allnames);
    n = 1;
    for k = 1:length(names)
        id = strmatch(names{k},allnames);
        name = names{k};
        ex = GetExptNumber(name);
        fprintf('Loading %s\n',name);
            if checkmode == 0 && sum([X(id).probes]) > 0
                load(name);
            end
        for j = id'
            ex = X(j).exptid;
            for k = 1:length(X(j).probes)
                cmd = sprintf(['AllVPcs(FullV, ''tchan'',' num2str(X(j).probes) ',' sprintf('%s, ', X(j).args{:}) ')']);
                fprintf([cmd '\n']);
                if checkmode == 0
                    if length(logfile)
                        logfid = fopen('logfile','a');
                        if logfid > 0
                            fprintf(logfid,'%s %s\n',datestr(now),cmd);
                        end
                    end
                    res{n} = AllVPcs(FullV, 'tchan', X(j).probes(k), X(j).args{:}, addargs{:});
                    fclose all;
                    res{n}.name = name;
                    res{n}.args = X(j).args;
                    res{n}.probes = X(j).probes(k);
                else
                    res.cmds{n} =  cmd;
                    if ex > 0
                    res.map(ex, X(j).probes(k)) = 1;
                    end
                end
                n = n+1;
            end
        end
        if checkmode == 0 & isfield(res{n-1},'toplevel') %AllVPCs was called
            close(res{n-1}.toplevel);
            clear FullV;
        end
    end


 function X = CheckExptList(X)
             %if args starts by listing all probes, strip this - want to call
        %one at at time next
        if strcmp(X.args{1},'tchan') && length(X.args{2}) > 2
            X.args = {X.args{3:end}};
        end
        if strcmp(X.args{1},'reclassifyall') 
            X.args{1} = 'reclassify';
        end

     for j = 1:size(X.cls)
         exptno = [];
         for k = 1:size(X.cls{j})
             exptno(k) = X.cls{j}{k}.cluster.exptno;
         end
         X.exptlist(j) = median(exptno);
     end
            
function result = CheckReclassify(X, checktype, varargin)
    
    exptno = 0;
    j = 1;
    while j <= length(varargin)
        if strncmp(varargin{j},'exptno',6)
            j = j+1;
            exptno = varargin{j};
        end
        j = j+1;
    end
    if isfield(X,'cls')
    [exps, exi] = sort(X.exptlist);
        acts = [];
        for j = 1:length(X.cls)
            e = X.exptlist(j);
            e = find(exps == X.exptlist(j));
            for k = 1:length(X.cls{j})
                if ~isfield(X.cls{j}{k},'probes')
                    X.cls{j}{k}.probes = k;
                end
                if ~isfield(X.cls{j}{k},'exptno')
                    X.cls{j}{k}.exptno = j;
                end
                X.cls{j}{k}.exptid = e;
                if ~isfield(X.cls{j}{k},'args')
                    X.cls{j}{k}.args = X.args;
                end
            end
            res = CheckReclassify(X.cls{j},checktype);
            result.actimage(e,:) = res.actimage(end,:);
            result.dateimage(e,:) = res.ctime;
            result.auto(e,:) = res.auto;
            result.xclim(e,:) = res.xcl;
            result.name{e} = X.name{j};
            if isfield(res,'readmethod')
                result.readmethod(e,:) = res.readmethod;
            end
            acts = [acts res.acts];
        end
        result.acts = acts;
        if isfield(X,'prefix')
            result.prefix = X.prefix;
        end
        return;
    end
    if iscell(X)
        res = CellToStruct(X);
    else
        f = {'probes' 'exptno'};
        for j = 1:length(X)
            for k = 1:length(f);
                if isfield(X,f{k});
                    res(j).(f{k}) = X(j).(f{k});
                end
            end
        end
    end
    for j = 1:length(X) 
        if iscell(X)
            res(j).result = CheckResult(X{j}, checktype);
            if isfield(X{j},'savetime')
            res(j).ctime = X{j}.cluster.savetime(1);
            else
            res(j).ctime = X{j}.cluster.ctime;
            end
            res(j).auto = X{j}.cluster.auto;
            if isfield(X{j}.cluster,'excludetrialids')  && length(X{j}.cluster.excludetrialids) > 1
                res(j).xcl = 1;
            else
                res(j).xcl = 0;
            end
            if isfield(X{j}.cluster,'exptreadmethod')
                res(j).readmethod = X{j}.cluster.exptreadmethod;
            end
            res(j).auto = X{j}.cluster.auto;
        else
            res(j).result = CheckResult(X(j),checktype);
        end
        
        if isfield(res,'args')
        if res(j).result == 1 %autocut
            res(j).args = {res(j).args{:} 'autocut' 'savespikes'};
        elseif res(j).result == 2
            res(j).args = {res(j).args{:} 'savespikes'};
        elseif res(j).result == 3
            res(j).args = {res(j).args{:} 'savespikes'};
        elseif res(j).result == 4
            res(j).args = {res(j).args{:} 'savespikes'};
        else 
            res(j).args = {};
        end
        end
        if isfield(res(j),'name')
            res(j).exptno= GetExptNumber(res(j).name);
        end
        actimage(res(j).exptno,res(j).probes) = res(j).result+1;
    end
    result.acts = res;
    result.actimage = actimage;
    result.ctime = [res.ctime];
    result.auto = [res.auto];
    result.xcl = [res.xcl];
    if isfield(res,'readmethod')
    result.readmethod = [res.readmethod];
    end
    
function res = CheckResult(a, type)
    if a.auto == 2 %old method for recording plotcluster cuts.
        a.auto = 2;
    end
      if strcmp(type,'reclassify')
          if a.err >0 && a.auto == 1
              res = 1;
          elseif a.err > 0
              fprintf('E%dP%d  err, not auto\n',a.cluster.exptno,a.cluster.probe(1)) 
              res = 0;
          elseif abs(a.xcorr(1)) > 0.95 && abs(a.xcorr(2)) > 0.95 && a.matchcounts(1)./a.matchcounts(2) < 0.15 && a.overlapn./a.matchcounts(2) > 0.5
              res = 4-a.auto; %4 or 3 or 2
          elseif a.auto ==1
              res = 1;
          elseif a.matchcounts(1) ==0 && a.matchcounts(3) == 0
              res = 5;
          else
              fprintf('E%dP%d xc %.2f %.2f overlap %d(%s)\n',...
                  a.cluster.exptno,a.cluster.probe(1),a.xcorr(1),a.xcorr(2),...
                  a.overlapn,sprintf(' %d',a.matchcounts));
              res = 0;
          end
      elseif strcmp(type,'baddim')
          if a.err ==1 && a.auto == 1
              res = 1;
          elseif a.err ==1 && a.auto == 0
              res = -1;
          else
              res = 0;
          end
      end
        