function [AllClusters, AllFullVData, details, AllClusterDetails] = LoadCluster(dirname, eid, varargin)
%LoadCluster(dirname, expts, ...)
%Load Cluster info for an expt or list of expts
%if length(expts) = 1 returns a 1 x nprobes cell matrix, else 
%a 1 x length(expts) matrix each containing 1 x nprobes cells
%LoadCluster(dirname, expts, 'array') keeps a cell arrray of cell arrays
%even if length(expts) ==1
%
%          ...,'nodetails') loads only the ClusterTimes file (quick)
%
%LoadCluster Combines info from Clusters and ClusterDetails files
%so that clst, t and Evec fields are in clusters
%LoadCluster(dirname, expts, 'gextxy') also inlcudes xy
%LoadCluster(dirname, expts, 'rawxy') also inlcudes xy with any roation
%removed (PlotClusters needs it this way
%LoadCluster(dirname, expts, 'alltimes') 
%replaces Clusters{}.times with ClusterDetails{}.t, saving memory
%(Also needed by PlotClusters)

AllClusters = {};
AllClusterDetails = {};
AllFullVData = {};
getauto = 0;
f = {'Evec' 'clst' 't'}; %fields to copy from Details
extraf = {'autofits'}; %fields to copy if absent
fixerrs = 0;
verbose = 0;
savespace = 0;

if ~isdir(dirname) && exist(dirname) %given filename
    name = dirname;
    dirname = fileparts(name);
    if nargin > 1
        varargin = {eid varargin{:}};
    end
    eid = GetExptNumber(name);
    if strfind(name,'AutoClusterTimes')
        getauto = 2;
    end
end



rawxy = 0;
alltimes = 0;
keeparray = 0;
getdetails = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'autoonly',8)
        getauto = 2;
    elseif strncmpi(varargin{j},'auto',4)
        getauto = 1;
    elseif strncmpi(varargin{j},'array',5)
        keeparray = 1;
    elseif strncmpi(varargin{j},'condense',5)
        f = {f{:} 'xy' 'triggerV' 'autofits'};
        savespace = 1;
    elseif strncmpi(varargin{j},'plain',5)
        plain = 1;
    elseif strncmpi(varargin{j},'fixerrs',5)
        fixerrs = 2;
    elseif strncmpi(varargin{j},'noauto',5)
        getauto = -1;
    elseif strncmpi(varargin{j},'nofix',5)
        fixerrs = 0;
    elseif strncmpi(varargin{j},'needdetails',7)
        getdetails =2;
    elseif strncmpi(varargin{j},'nodetails',5)
        getdetails =0;
    elseif strncmpi(varargin{j},'quick',5)
        getdetails =0;
    elseif strncmpi(varargin{j},'getxy',5)
        f = {f{:} 'xy' 'triggerV' 'autofits'};
    elseif strncmpi(varargin{j},'rawxy',5)
        rawxy = 1;
        f = {f{:} 'xy' 'triggerV'};
    elseif strncmpi(varargin{j},'alltimes',5)
        alltimes = 1;
    elseif strncmpi(varargin{j},'savespace',5)
        savespace = 1;
    elseif strncmpi(varargin{j},'verbose',5)
        verbose =1;
    end
    j = j+1;
end
details = {};

dirname = CheckNameBug(dirname);

if ~isdir(dirname)
    cprintf('red','%s is not a Diretcory\n',dirname);
    return;
end

loadname = '';
for j = 1:length(eid)
    ts = now;
    e = floor(eid(j));
    exptnos(j) = e;
    if round(rem(eid(j),e).*10) == 1
        xs = 'a';
    else
        xs = '';
    end
    name = [dirname '/Expt' num2str(e) xs 'ClusterTimes.mat'];
    dname = [dirname '/Expt' num2str(e) xs 'ClusterTimesDetails.mat'];
    daname = [dirname '/Expt' num2str(e) xs 'AutoClusterTimesDetails.mat'];
    aname = [dirname '/Expt' num2str(e) xs 'AutoClusterTimes.mat'];
    if getauto > 0
        dname = daname;
        name = aname;
    end
        
    needf = 0;
    if exist(aname,'file') && getauto >= 0
        names{j} = aname;
        tic;
        try
        load(aname);        
        details{j}.loaddur = toc;
        catch
            cprintf('red','Error Loading %s\n',aname);
            Clusters = {};
        end
        d = dir(aname);
        details{j}.loadbytes = d.bytes;
        
        AutoClusters = Clusters;
        for p = 1:length(AutoClusters)
            AutoClusters{p}.auto = 1;
            AutoClusters{p}.loadname = aname;
            if isfield(AutoClusters{p},'eckercluster') && isfield(AutoClusters{p}.eckercluster,'fit')
                AutoClusters{p}.eckercluster.fit = rmfields(AutoClusters{p}.eckercluster.fit,'Y');
            end
            for k = 1:length(f) 
               if ~isfield(AutoClusters{p},f{k})
                needf = needf+1;
               end
            end
        end
        if getdetails == 2 %force loading of details
            needf = 1;
        end
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
        details{j}.loadname = aname;
    else 
        names{j} = name;
        AutoClusters = {};
    end
    if getauto == 2
        AllClusters{j} = AutoClusters;
        fprintf('Loading %s at %s\n',name,datestr(now,'HH:MM:SS'));
    elseif exist(name,'file')
        tic;
        load(name);
        names{j} = name;
        details{j}.loaddur(1) = toc;
        d = dir(name);
        details{j}.loadbytes(1) = d.bytes;
        
        details{j}.loadname = name;
        AllClusters{j} = Clusters;
        for p = 1:length(Clusters)
            for k = 1:length(f)
                if ~isfield(Clusters{p},f{k})
                    needf = needf+1;
                end
            end
        end
        loadname = name;
        AllClusters{j}{1}.loadname = name;
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
        for k = 1:length(AutoClusters)
            if k > length(Clusters) || isempty(Clusters{k})
                AllClusters{j}{k} = AutoClusters{k};
            end
        end
    elseif ~isempty(AutoClusters)
        AllClusters{j} = AutoClusters;
        details{j}.err = sprintf('Can''t read %s\n',name);
        fprintf('Can''t read %s, but loaded %s (%.2f sec)\n',name,aname,mytoc(ts));
    else
        fprintf('Can''t read %s or %s (%.2f sec)\n',name,aname,mytoc(ts));
        AllClusters{j} = {};
        AllFullVData{j} = [];
        details{j}.err = sprintf('Can''t read %s or %s (%.3f sec)\n',name,aname,mytoc(ts));
    end
    details{j}.loadtime = now;
    details{j}.exptno = e;
    if getdetails && needf && isfield(details{j},'loaddur') %did load ClusterTimes
        [ClusterDetails, CD] = LoadClusterDetails(name);
        details{j}.loaddur = cat(2,details{j}.loaddur, CD.loaddur);
        details{j}.loadbytes = cat(2,details{j}.loadbytes, CD.loadbytes);
        for k = 1:length(ClusterDetails)
            details{j}.probe = k;
            C = FixCluster(AllClusters{j}{k});
            CD = ClusterDetails{k};
            ClusterDetails{k}.exptno = C.exptno;
            if isfield(ClusterDetails{k},'clst')
                if length(ClusterDetails{k}.clst) ~= AllClusters{j}{k}.nspks
                    details{j} = AddError(details{j},'Cluster %d E%d Details Clst (%d) does not match  Cluster nspks (%d)!!\n',k,e,length(ClusterDetails{k}.clst),AllClusters{j}{k}.nspks);
                    AllClusters{j}{k}.needed = 1;
                elseif sum(ClusterDetails{k}.clst ==2) ~= length(AllClusters{j}{k}.times)
                    cid = find(CD.clst ==2);
                    a = length(cid);
                    b = length(AllClusters{j}{k}.times);
                    if isfield(C,'clusterprog')
                        pstr = C.clusterprog;
                    else
                        pstr = 'AllVPcs';
                    end
                    useclst = 0;
                    ncut = AllClusters{j}{k}.ncut;
                    useclsts = [];
                    if C.ncut == length(C.times)
                        useclsts(1) = 2;
                    else
                        useclsts(1) = 0;
                    end
                        
                    for c = 1:length(C.next)
                        if isfield(C.next{c},'ctime') && C.next{c}.ctime > max(C.savetime)
                            if isfield(C.next{c},'clusterprog')
                                cstr = C.next{c}.clusterprog;
                            else
                                cstr = 'Unknown';
                            end
                            if strncmp(cstr,'PlotClusters',10) & C.next{c}.ncut == length(C.next{c}.times)
                                fprintf('E%dP%dCluster%d last saved by %s. Setting Details.clst from times',C.exptno,k,c+1,cstr);
                                [ta,mid,tb,cmid] = MatchTimes(CD.t,C.times,0.0005);
                                CD.clst(mid) = c+2;
                                useclst = 2;
                                useclsts(c) = 2;
                            elseif ncut ~= b
                                useclst =1;
                                fprintf('E%dP%dCluster%d (by %s) newer than Save!! Using clst for times',C.exptno,k,c+1,cstr);
                                AllClusters{j}{k}.fixmode = 'useclst';
                                AllClusters{j}{k}.oldtimes = C.times;
                                AllClusters{j}{k}.times = CD.t(CD.clst ==2);
                                AllClusters{j}{k}.ncut = sum(CD.clst ==2);
                            else
                                fprintf('E%dP%dCluster%d (by %s) newer than Save!!',C.exptno,k,c+1,cstr);
                            end
                        end
                    end
                    if sum(useclsts ==2)
                        fprintf('Fixing clst with Cluster.times\n');
                        AllClusters{j}{k}.fixmode = 'usetimes';
                        AllClusters{j}{k}.oldcounts = Counts(CD.clst);
                        CD.clst(1:end) = 1;
                        for c =0:length(C.next)
                            if c == 0
                                [ta,mid,tb,cmid] = MatchTimes(CD.t,C.times,0.0005);
                            elseif isfield(C.next{c},'times')
                                [ta,mid,tb,cmid] = MatchTimes(CD.t,C.next{c}.times,0.0005);
                            else
                                mid = [];
                            end
                            CD.clst(mid) = c+2;
                        end
                    elseif useclst == 0
                        fprintf('\n');
                        details{j} = AddError(details{j},'Cluster %d E%d (by%s on %s) Details Mismatch: Clst 1 is %d, length(Cluster.times) is %d (ncut is%d) Details saved %s!!\n',...
                            k,e,pstr,datestr(C.savetime(1)),a,b,AllClusters{j}{k}.ncut,datestr(CD.savetime));
                        AllClusters{j}{k}.needed = 1;
                        if C.savetime(1) > C.savetime(2) && strncmp(pstr,'PlotClusters',8)
                            cprintf('blue', 'Was quantifyall run without ''savespikes'' ?\n');
                        end
                    end
                    if fixerrs
                        [a,mid,b,cmid] = MatchTimes(CD.t(cid),C.times,0.0005);
                        [a,id,b,c] = MatchTimes(C.times,CD.t(cid),0.0005);
                        if fixerrs == 2 %plot results
                            GetFigure('Missing events');
                            hold off;
                            plot(C.times(b),CD.t(cid(a)),'.'); %matches
                            hold on;
                            plot(C.times(c),CD.t(cid(id)),'r.'); %clst not in times
                            plot(C.times(mid),CD.t(cid(cmid)),'g.'); %times not in clst, nearest clst
                            xlabel('ClusterTimes')
                            ylabel('clst==2');
                            legend('matches','clst not in times', 'times not in clst')
                            [a,~,b,c] = MatchTimes(C.times(mid),CD.t,0.0002); %times in clst but not classified
                            plot(C.times(mid(b)),CD.t(a),'+');
                            if length(b) < length(mid)
                                fprintf('%d times not in clst\n',length(mid)-length(b));
                            else
                                fprintf('All times in clst,%d not cluster1\n',length(mid));
                            end
                        end
                        AllClusters{j}{k}.times = sort([C.times CD.t(cid(id))]);
                    end
                end
            elseif isfield(CD,'t') && isfield(C,'times')
               ClusterDetails{k}.clst = ones(length(CD.t),1);
               id = find(ismember(CD.t,C.times));
               ClusterDetails{k}.clst(id) = 2;
            end
            if isfield(ClusterDetails{k},'ctime') && ...
                    isfield(AllClusters{j}{k},'savetime') && ...
                    ClusterDetails{k}.ctime > AllClusters{j}{k}.savetime(end) + 0.0001; %10 sec diff
                fprintf('Cluster %d Details newer than Cluster!!\n',k);
                AllClusters{j}{k}.needed = 1;
            end
            for n = 1:length(f)
                if isfield(ClusterDetails{k},f{n})
                    AllClusters{j}{k}.(f{n}) = ClusterDetails{k}.(f{n});
                end
            end
            for n = 1:length(extraf)
                if isfield(ClusterDetails{k},extraf{n}) && ~isfield(AllClusters{j}{k},extraf{n})
                    AllClusters{j}{k}.(extraf{n}) = ClusterDetails{k}.(extraf{n});
                end
            end
            if isfield(ClusterDetails{k},'next')
                for c = 1:length(ClusterDetails{k}.next)
                    %only load ClusterDetails if Cluster is still defined in next{c}
                    if c <= length(AllClusters{j}{k}.next) && isfield(AllClusters{j}{k}.next{c},'space')
                        for n = 1:length(f)
                            if isfield(ClusterDetails{k}.next{c},f{n})
                                AllClusters{j}{k}.next{c}.(f{n}) = ClusterDetails{k}.next{c}.(f{n});
                            end
                        end
                    end
                end
            end
%if no xy in next, AND next{c} i snot emptyshould mean it used the same as Cluster 1            
            for c = 1:length(AllClusters{j}{k}.next) 
                if ~isfield(ClusterDetails{k},'next') || c > length(ClusterDetails{k}.next) || ~isfield(ClusterDetails{k}.next{c},'xy')
                    if isfield(ClusterDetails{k},'xy') && isfield(AllClusters{j}{k},'xy') ... 
                            && ~isempty(ClusterDetails{k}.xy) && ~isempty(AllClusters{j}{k}.next{c})
                        AllClusters{j}{k}.next{c}.xy = ClusterDetails{k}.xy;
                    end
                end
            end
            %some old files saved as row, not column
            if isfield(AllClusters{j},'clst') && size(AllClusters{j}(k).clst,1) == 1
                AllClusters{j}(k).clst = AllClusters{j}(k).clst';
            end
            if rawxy && isfield(ClusterDetails{k},'xy')
                xy = ClusterDetails{k}.xy;
                C = AllClusters{j}{k};
                if C.shape == 0
                    AllClusters{j}{k}.xy = xy;
                else
                    AllClusters{j}{k}.xy = xyrotate(xy(:,1),xy(:,2),-C.angle);
                end
            end
            if alltimes && isfield(AllClusters{j}{k},'t')
                AllClusters{j}{k}.times = AllClusters{j}{k}.t;
                AllClusters{j}{k} = rmfield(AllClusters{j}{k},'t');
            end
        end
        AllClusterDetails{j} = ClusterDetails;
    else
        details{j}.loaddur(2:3) = 0;
    end
end

for j = 1:length(AllClusters)
    name = names{j};
    for k = 1:length(AllClusters{j})
        if ~isempty(loadname)
            AllClusters{j}{k}.loadname = loadname;
        end
%         C = FixCluster(AllClusters{j}{k});
        C = AllClusters{j}{k};
        if savespace
            if isfield(C,'eckercluster')
                C.eckercluster = rmfields(C.eckercluster,'fit');
            end
            if isfield(C,'Evec') && isfield(C.Evec,'Evec')
                C.Evec.Evec = C.Evec.Evec(1:6,:);
            end
            if isfield(C,'autofits')
                for c = 1:length(C.autofits)
                    C.autofits{c} = rmfields(C.autofits{c},'eckercluster');
                end
            end
            C.condensed = 1;
        end
        if isfield(C,'space')
            if verbose
                fprintf('Checking E%dP%d\n',j,k);
            end
            if ~isfield(C,'quick')
                C.quick = 0;
            end
            if ~isfield(C,'dropi')
                C.dropi = [0 0 0 0];
            end
            if ~isfield(C,'trigdt')
                C.trigdt = 0;
            end
            if ~isfield(C,'manual')
                C.manual = 0;
            end
            if ~isfield(C,'next')
                C.next = {};
            end
            if ~isfield(C,'clusterprog')
                C.clusterprog = '';
            end
            if ~isfield(C,'next')
                C.next = {};
            elseif isstruct(C.next)
                next = C.next;
                C = rmfield(C,'next');
                C.next{1} = next;
            end
            AllClusters{j}{k} = C;
            if isfield(C,'trigset')
                for c = 1:length(C.trigset{1}.next)
                    if c > length(C.next) || isempty(C.next{c}) || (isfield(C.next{c},'triggerset') &&C.next{c}.triggerset ==1)
                        if isfield(C.trigset{1}.next{c},'space')
                            C.next{c} = C.trigset{1}.next{c};
                            C.next{c}.trigset = 1;
                        end
                    end
                end
                AllClusters{j}{k} = rmfields(C, 'pcplot');
                if isfield(ClusterDetails{k},'trigset') && isfield(ClusterDetails{k}.trigset{1},'clst')
                    AllClusters{j}{k}.trigset{1}.clst = ClusterDetails{k}.trigset{1}.clst;
                    AllClusters{j}{k}.trigset{1}.t = ClusterDetails{k}.trigset{1}.t;
                    Cd = ClusterDetails{k}.trigset{1};
                    for c = 1:length(Cd.next)
                        if isfield(Cd.next{c},'xy') && (c > length(AllClusters{j}{k}.next) || ~isfield(AllClusters{j}{k}.next{c},'xy'))
                            AllClusters{j}{k}.next{c}.xy = Cd.next{c}.xy;
                        end
                    end                    
                end
            end
                if isfield(AllClusters{j}{k},'times') && diff(size(AllClusters{j}{k}.times)) < 0
                    if isempty(AllClusters{j}{k}.times)
                        cprintf('red','%s P%d times is empty\n',name,k);
                    else
                        cprintf('red','%s P%d times is a row vector\n',name,k);
                        AllClusters{j}{k}.times = AllClusters{j}{k}.times';
                    end
                end
                if isfield(AllClusters{j}{k},'t') && diff(size(AllClusters{j}{k}.t)) < 0
                    cprintf('red','%s P%d t is a row vector\n',name,k);
                    AllClusters{j}{k}.t = AllClusters{j}{k}.t';
                end
            end
            details{j}.exptno = eid;
    end
end

if length(AllClusters) == 1 && keeparray == 0
    AllClusters = FixCluster(AllClusters{1});
    if ~isempty(AllClusterDetails)
        AllClusterDetails = AllClusterDetails{1};
    end
    if length(AllFullVData)
    AllFullVData = AllFullVData{1};
    end
    details = details{1};
end


