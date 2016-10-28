function [DATA, finddetails] = FindCells(DATA, varargin )
%DATA = PC.FindCells(DATA) Attempts to fill in CellList automatically
%Adds field DATA.AutoCellList
%PC.FindCells(DATA, DATA.AutoCellList) re-does looking for duplicates, and
%PC.FindCells(DATA, 'rebuild') forces start from scratch
%PC.FindCells(DATA, 'rebuild', 'save') saves list and details at the end
%PC.FindCells(DATA, 'stage',n) forces start stage n if length(n) is 2,
%                                                      stops at stage n(2);
%PC.FindCells(DATA, 'stopstage',n) set stopping point
%
%    ...,'nofixes') stops changes being made to Clusters
%    ...,'applyfixes') changes that impove cell matches are applied Clusters
%                      this is the default
% to plot the results use PC.PlotCellList(DATA,'autolist','showfig') 
%PC.FindCells(DATA, 'set') Copies the CellList and modified Cluster data
%  so that the plot shows the final list
%PC.FindCells(DATA, 'set',n ) set Data to match stage n
%
%after [xDATA, details] = PC.FindCells
%PC.FindCells(DATA,xDATA, details) or
%PC.FindCells(figurehandle,xDATA, details) 
%updates DATA.autolist, and
%AutoClusters, so that can use DATA for plotting auto list 
%
%
%Stage 1: Initial list based on finding good cells and follwoing them
%Stage 2: Remove Duplicate
%Stage 3: Find Orphan squares/ delete isolated squares
%Stage 4: now empty
%Stage 5: merge vertically separted cells with same shape
%Stage 6: Old method for 5
%Stage 7: Delete Isolated squares
%Stage 8: Remove Duplicates
%Stage 9: Fix Boundary Changes
%Stage 10: Look for "Gaps" in a cell and check they cannot be fille
%Stage 11: Tidy up. Re-order cell numbers, swap clusters around so cell has
%same cluster# where poss, 
%
%PC.FindCells(DATA, 'findcell',[e p cl; eb 0 0]) call routine to 
%             find cluster in expt eb, tha matches [e p cl]. For testing
%PC.FindCells(DATA, 'comparecells',[a b])  Checks if a should merge with b
%PC.FindCells(DATA, 'comparecells',[a a])  does self consistency check for cell a
%
%PC.FindCells(DATA, 'checkshapes',stage)  Checks if a should merge with b
%          Checks each cell in CellLists(stage) for shape consistency
%
%Bruce Notes
%PC.FindCells(DATA,'stage',100) runs and applies celllist.FixDuplicates
%  'saveinit') Saves 'InitialClusters.mat' = before any fixing. 

%To do 
% look at the  dropcells list automatically to see if there is a
% consistent cell being dropped.  Not just MU of a real cell or just random
%pprM034 P19 E15-40 looks like a good example. Might be worth retriggering
%even....
%
% Candidates. Allow recalcuation of this iwth changing thresholds to see
% what is best.
%? Show just the best, and if its clicked, show the next one?
verbose = 1;
CellList = [];
celldup = [];
rebuild = 0;
listeid = 1;
stopstage = NaN;
recordstate = 100;
force.buildduplicates = 0;
force.addmu = 0;
force.stage = -1;
force.swapclusters = -1; %don't mess with this until everything else works
force.noswap = 0;
checkshapes = 0;
state.watchprogress = 1;
state.verbose = 1;
state.applyfixes = 1;
alignclusters = 0; %tries to make cluster number uniform for a cell

%These are the parameter that determine the initial search
iscrit = [3.1 2.8];
dropcrit = 1.8;
useauto =1;
manualfind = [];
comparecells = [];
stopcell = NaN;

D = [];
ncheck = 1;
if length(varargin) > 1  && isfield(varargin{1},'autolist') && isfield(varargin{2},'Clusters')
    if isfigure(DATA)
        F = DATA;
        DATA = get(F,'UserData');
    end
    if isfield(DATA,'SpaceTypeLabels')
    DATA.autolist = varargin{1}.autolist;
    setappdata(DATA.toplevel,'AutoClusters',varargin{2}.Clusters);
    SetData(DATA);
    DATA = PC.PlotCellList(DATA,'autolist','showfig');
    return;
    end
end
if isfield(DATA,'Q')  %Given an existing finddetails: review contents
    CheckDetails(DATA, varargin{:});
    return;
end

if isfigure(DATA)
    DATA = GetDataFromFig(DATA);
end

DATA.autolist.listtype = 'autolist';
Clusters = PC.GetValue(DATA,'Clusters','withfits');
DATA.tag = AddField(DATA.tag,'autocelllist', ['auto' DATA.tag.celllist]);

strargs = cell2cellstr(varargin);
Expts = getappdata(DATA.toplevel,'Expts');

if ~isempty(varargin)
    if sum(strcmp(varargin{1},{'set'}))
        if force.stage < 0
            maxstage = 12;
        else
            maxstage = force.stage;
        end
        finddetails = getappdata(DATA.toplevel,'AutoCellList');
        if length(varargin) > 1 && isnumeric(varargin{2}) && varargin{2} > 0
            stage = varargin{2};
        else
            stage = length(finddetails.CellLists);            
        end
        DATA.autolist.CellList = finddetails.CellLists{stage};
        if ~isfield(finddetails,'clusterswaps')
            finddetails.clusterswaps = [];
        end
        j = 1;
        while isempty(DATA.autolist.CellList) && j < stage
            fprintf('Empty CellList at Stage %d. Trying %d\n',stage,stage-j);
            DATA.autolist.CellList = finddetails.CellLists{stage-j};
            j = j+1;
        end
            
        finddetails.CellList = finddetails.CellLists{stage};        
%reverse any swaps first just in case.        
        if isfield(finddetails,'ChangedClusters')
            for j = 1:length(finddetails.ChangedClusters)
                if j <= maxstage && ~isempty(finddetails.ChangedClusters{j})
                    Clusters = clust.Merge(Clusters,finddetails.ChangedClusters{j});
                end
            end
        end
        Clusters = clust.Swap(Clusters,'reverse');
        if stage ==1
            if isappdata(DATA.toplevel,'DriftFixClusters')
                Clusters = getappdata(DATA.toplevel,'DriftFixClusters')
            else
                [X, Clusters] = PC.CheckClusters(Clusters,'autofits');
                [Clusters, cdetails] = clust.Check(Clusters,'drift','apply', DATA.ArrayConfig,state,Expts);
%                setappdata(DATA.toplevel,'DriftFixClusters',Clusters);
                DATA.autolist.driftfix = 1;
            end
        elseif stage < 11
            setappdata(DATA.toplevel,'AutoClusters',Clusters);
        else
            if isfield(finddetails,'muCellList')
                DATA.autolist.muCellList = finddetails.muCellList;
            end
            Clusters = clust.Swap(Clusters,finddetails.clusterswaps);
            setappdata(DATA.toplevel,'AutoClusters',Clusters);
        end
%rebuild quality measures to refelct changes in Clusters        
       if (finddetails.stage >= 11 && stage < 11) || (finddetails.stage <11 && stage >= 11)
           [Q, D, I, cellA] = PC.MakeVarMatrix(DATA,Clusters,'cellquality','dropi','defaultisolation','probeamplitude');
           finddetails.Q = Q; finddetails.D = D; finddetails.I = I; finddetails.A = cellA;
           DATA.autolist.dropi = D;
       end
       DATA.autolist = CopyFields(DATA.autolist,finddetails, {'corrproberange' 'corrth'});
       finddetails.stage = stage;
        setappdata(DATA.toplevel,'AutoCellList',finddetails);
        DATA.autolist.listtype = 'autolist';
        SetData(DATA);
        if strcmp(varargin{j},'set')
%tag.autocelllist will be set            
            DATA = PC.PlotCellList(DATA,'autolist','showfig');
            return;
        end
    end
    
    if strcmp(varargin{1},'reset')
        finddetails = getappdata(DATA.toplevel,'AutoCellList');
        DATA.autolist.CellList = finddetails.CellLists{10};
        DATA.autolist.dropi = finddetails.D;
        setappdata(DATA.toplevel,'AutoClusters',finddetails.StartClusters);
        SetData(DATA);
        PC.PlotCellList(DATA,'autolist','showfig');
        return;
    end
    if strncmp(varargin{1},'save',4)
        SaveAutoList(DATA);
        if strncmp(varargin{1},'saveall',7);
            savename = sprintf('%s/FixedClusters.mat',DATA.autodatadir);
            save(savename, 'Clusters');
        end
        return;
    end
    if strcmp(varargin{1},'load')
%Load should be done when Origincal Clusters are in place. And Should then
%run FindCells('stage', 10) to make swapped clusters
        savename = sprintf('%s/AutoCellList.mat',DATA.autodatadir);
        d = dir(savename);
        fprintf('Loading %s (%s)\n',savename,d.date);
        X = my.load(savename);
        finddetails = X.finddetails;
        if ~isfield(finddetails,'stage')
            finddetails.stage = 12; %assume finished
        end
        finddetails = AddField(finddetails,'corrproberange',2);
        setappdata(DATA.toplevel,'AutoCellList',finddetails);
        
        [a,b] = cellstrcmp('set',varargin);
        if a
            a = find(b >0);
            if length(varargin) > a && isnumeric(varargin{a+1})
                force.stage = varargin{a+1};
            end
        end
        X = my.load([DATA.autodatadir '/InitialFixedClusters.mat'],'safe');
        if isfield(X,'ChangedClusters')
            AClusters = Clusters;
            Clusters = clust.Merge(Clusters, X.ChangedClusters);
            finddetails.ReverseFix = clust.Compare(Clusters, AClusters,'diffonly','noplot');
            clear AClusters;
            setappdata(DATA.toplevel,'AutoClusters',Clusters);
        end
        
        savename = sprintf('%s/FixedClusters.mat',DATA.autodatadir);
        runt = max(finddetails.t(1:5));
        d = dir(savename);
        if ~isempty(d) && d.datenum > runt-1 && (force.stage <0 || force.stage > 5)%saved less than a day before run
            X = my.load(savename);
        end
        if sum(strcmp('set',strargs))
            DATA = PC.FindCells(DATA.toplevel,'set',force.stage);
        end
        return;
    end
    if strcmp(varargin{1},'findcell')
    end
end


    
corrth = [0.97 0.94];
minexpts = 3;  %cells only present for fewer expts are removed
%find expts we want to ingore. skipeid is a list of rows, not exptnos
ed = prctile(DATA.electrodedepth,50);
DATA.autolist.eid = find(abs(DATA.electrodedepth(DATA.exptlist) - ed) <= 0.1);
DATA.autolist.skipeid = setdiff(1:length(DATA.exptid),DATA.autolist.eid);

fixduplicates = 0;
if isappdata(DATA.toplevel,'AllDuplicates')
    buildduplicates = 0;
else
    buildduplicates = 1;
end

j = 1;
while j <= length(varargin)
    if isnumeric(varargin{j}) && ndims(varargin{j}) ==3
        CellList = varargin{j};
    elseif strcmp(varargin{j},'applyfixes')
        state.applyfixes = 1;
    elseif strcmp(varargin{j},'savefixes')
        ApplyFixes(DATA);
        return;
    elseif strcmp(varargin{j},'buildduplicates')
        force.buildduplicates  = 1;
    elseif strncmp(varargin{j},'checkshapes',7)
        j = j+1;
        checkshapes = varargin{j};
    elseif strncmp(varargin{j},'compare',7)
        j = j+1;
        comparecells = varargin{j};
    elseif strcmp(varargin{j},'expts')
        j = j+1;
        DATA.autolist.eid = varargin{j};
    elseif strcmp(varargin{j},'nofixes')
        state.applyfixes = 0;
    elseif strcmp(varargin{j},'nomu')
        force.addmu = -1;
    elseif strcmp(varargin{j},'nowatch')
        state.watchprogress = 0;
    elseif strcmp(varargin{j},'noswap')
        force.swapclusters = -1;
    elseif strcmp(varargin{j},'stage')
        j = j+1;
        force.stage = varargin{j}(1);
        if length(varargin{j}) > 1
            stopstage = varargin{j}(2);
        end
    elseif strcmp(varargin{j},'stopstage')
        j = 1+j;
        stopstage = varargin{j};
    elseif strcmp(varargin{j},'useauto')
        useauto = 1;
    elseif strncmp(varargin{j},'finddup',7)
        fixduplicates = 2;
    elseif strcmp(varargin{j},'findcell')
        j = j+1;
        manualfind = varargin{j};
        listeid = 0;
    elseif strcmp(varargin{j},'stopcell')
        j = j+1;
        stopcell = varargin{j};
    elseif strncmp(varargin{j},'rebuild',7)
        rebuild = 1;
        if strcmp(varargin{j},'rebuildall')
            rebuild = 2;
            if isfield(DATA.autolist,'driftfix') && DATA.autolist.driftfix > 0;
            tmpname = [DATA.autodatadir '/InitialClusters.mat'];
            if exist(tmpname)
                tic;
                load(tmpname);
                toc
            end
            DATA.autolist.driftfix = 0;
            end
        end
        buildduplicates = 1;
        findduplicates = 1;
    elseif ischar(varargin{j})
        strargs = {strargs{:} varargin{j}};
    end
    j = j+1;
end

if force.stage

end

if ~isempty(DATA.autolist.skipeid) && listeid
    cprintf('red','Skpping Expts with electrode depth different from %.3f:\n',ed);
    for j = 1:length(DATA.autolist.skipeid)
        e = DATA.autolist.skipeid(j);
        fprintf('%d(%s):%.2f\n',e,DATA.expnames{e},DATA.electrodedepth(e));
    end
end

finddetails.toplevel = DATA.toplevel;
if fixduplicates > 1
    FindDoubleTriggers(Clusters);
end
checkcell(1).eid = 0;

if fixduplicates ==2 %test of shape comparison over all expts
   CompareCellShape(DATA,Clusters,DATA.autolist.CellList,3,6);
end


if force.buildduplicates || ~isappdata(DATA.toplevel,'AllDuplicates')
    buildduplicates =1;
    fixduplicates = 1;
end
    

    
finddetails.t(1) = now;

if isempty(CellList) && isappdata(DATA.toplevel,'AutoCellList') && rebuild == 0
    x = getappdata(DATA.toplevel,'AutoCellList');
    if isfield(x,'CellList')
        CellList = x.CellList;
    elseif isfield(x,'InitialList')
        CellList = x.InitialList;
    end
    fprintf('Using AuotCellList stored in appdata\n');
    fixduplicates = 1;
    finddetails = x;
    finddetails.toplevel = DATA.toplevel;
    if isfield(finddetails,'checkcell')
        ncheck = length(finddetails.checkcell);
        checkcell = finddetails.checkcell;
    else
        ncheck = 1;
    end
    if force.stage >= 0
        finddetails.stage = force.stage;
    end
    if ~isfield(finddetails,'stage')
        finddetails.stage = 1;
    elseif finddetails.stage > 11 %extra manual steps/tests
        CellList = finddetails.CellLists{11};        
    elseif finddetails.stage == 0 %start over
        CellList = [];
    else
        CellList = finddetails.CellLists{finddetails.stage};        
        if isempty(CellList)
            CellList = finddetails.CellLists{finddetails.stage-1};
        end
    end
    if finddetails.stage > 9 && isappdata(DATA.toplevel,'NewClusters')
        Clusters = getappdata(DATA.toplevel,'NewClusters');
    end
%make this match the selected autolist in case plot stuff below    
    DATA.autolist.CellList = CellList;
    SetData(DATA);  
    if isfield(finddetails,'A')
        D = finddetails.D;
        I = finddetails.I;
        Q = finddetails.Q;
        cellA = finddetails.A;
    end
end

if ~isfield(finddetails,'ArrayConfig')
    finddetails.ArrayConfig = [];
end
if isfield(DATA,'ArrayConfig') && isempty(finddetails.ArrayConfig)
    finddetails.ArrayConfig = DATA.ArrayConfig;
    if isfield(DATA.ArrayConfig,'X') && length(unique(DATA.ArrayConfig.X)) == 2
        finddetails.corrproberange = 2; %up to 8 probes
    else
        finddetails.corrproberange = 5;
    end
end
finddetails.corrth = corrth;
DATA.autolist= CopyFields(DATA.autolist, finddetails,{'corrproberange' 'corrth'});

savetime = PC.MakeVarMatrix(DATA,Clusters,'2D','savetime');
if isfield(finddetails, 'savetime')
    newc = find(savetime > finddetails.savetime);
    if ~isempty(newc)
        fprint('%d Clusters changed since last build\n',length(newc));
    end
end

fixlist = [];
if isempty(D)
    LogState(DATA,finddetails);
    tmpname = [DATA.autodatadir '/InitialClusters.mat'];
    if ~exist(tmpname)  || cellstrcmp('saveinit',varargin);
        save(tmpname,'Clusters');
    end
%first find cases where drift has caused splitting into two clusters, and
%fix these
    iname = [DATA.autodatadir '/InitialFixedClusters.mat'];
    if ~isfield(DATA.autolist,'driftfix') || DATA.autolist.driftfix == 0
        if rebuild == 2 || ~exist(iname)%really starting over
        AClusters = Clusters;
        fprintf('Checking Auto Clusters\n');
        [X, Clusters] = PC.CheckClusters(Clusters,'autofits');
        fprintf('Applying Fixes\n');        
        [Clusters, cdetails] = clust.Check(Clusters,'drift','apply', DATA.ArrayConfig,state,Expts,'noninteractive');
        sz = size(cdetails.fixfit);
        fixlist = reshape(shiftdim(cdetails.fixfit,2),sz(1)*sz(3),sz(2));
        id = find(sum(fixlist')>0);
        fixlist = fixlist(id,:);
        fixlist(:,size(fixlist,2)+1) = 1;  %who did fix
        fixlist(:,size(fixlist,2)+1) = now;  %who did fix

        finddetails.shakelist = find(cdetails.isshake>0);

        sz = size(cdetails.combinecls);
        ccls = reshape(cdetails.combinecls,sz(1)*sz(2),sz(3));
        finddetails.driftcombine = ccls(id,1:2);
        DATA.autolist.initialfix = fixlist;
        DATA.autolist.shakelist = finddetails.shakelist;
        [ChangedClusters, details] = clust.Compare(AClusters,Clusters,'diffonly','noplot');
        clear AClusters;
        setappdata(DATA.toplevel,'AutoClusters',Clusters);
%        setappdata(DATA.toplevel,'DriftFixClusters',Clusters);
        save(iname,'ChangedClusters');
        else
            load(iname);
            Clusters = clust.Merge(Clusters,ChangedClusters);
            setappdata(DATA.toplevel,'AutoClusters',Clusters);
        end
        DATA.autolist.driftfix = 1;
        SetData(DATA);
    elseif isfield(DATA.autolist,'initialfix')
        if isappdata(1,'DriftFixClusters')
            Clusters = getappdata(DATA.toplevel,'DriftFixClusters');
        elseif exist(iname)
            tic;
            fprintf('Loading %s...',iname)
            X = load(iname);
            if isfield(X,'ChangedClusters')
                Clusters = clust.Merge(Clusters, X.ChangedClusters);
            end
            toc;
        end
        fixlist = DATA.autolist.initialfix;
        if isfield(DATA.autolist,'shakelist')
            finddetails.shakelist = DATA.autolist.shakelist;
        else
            finddetails.shakelist = [];
        end
    end
    [Q, D, I, cellA] = PC.MakeVarMatrix(DATA,Clusters,'cellquality','dropi','defaultisolation','probeamplitude');
    finddetails.D = D;
    finddetails.I = I;
    finddetails.Q = Q;
    finddetails.A = cellA;
    finddetails.corrth = corrth;
    [~, nccount] = clust.Check(Clusters,'ncluster')
    if isfield(nccount,'nc')
        finddetails.nc = nccount.nc;
    end

elseif finddetails.stage > 1
    if isappdata(DATA.toplevel,'DriftFixClusters')
%        Clusters = getappdata(DATA.toplevel,'DriftFixClusters');
    end
end
fprintf('Finished Building Data after %.2f sec\n',mytoc(finddetails.t(1)));
finddetails.verbose = verbose;
if ~isfield(finddetails,'stage')
    finddetails.stage = 0;
end
    
if finddetails.stage >= 100  %special case. Try fixing cases where cell and duplicate appear on a single probe
    C = getappdata(DATA.toplevel,'Clusterbackup')
    if length(C) == length(Clusters)
        Clusters = C; %start over
    end
    FixList(DATA,  Clusters, CellList, finddetails.stage);
    return;
end
if rebuild
    if isappdata(DATA.toplevel,'SelfDetails')
        rmappdata(DATA.toplevel,'SelfDetails');
    end
    LogState(DATA,finddetails);
end
if checkshapes
    x = CheckSelfShapes(DATA, finddetails.CellLists{checkshapes},Clusters, finddetails);
    PlotComparison(DATA,x,'fullwait');
    return;
end

if ~isempty(comparecells)
    x = PC.CompareCells(DATA, CellList,Clusters,comparecells,'plot');
    return;
end

if ~isempty(manualfind)
    e = manualfind(2,1);
     [C, details] = FindCellForExpt(DATA, e, Clusters, Q,D,I, manualfind(1,:),state, finddetails,'corrth',corrth,'wide',DATA,'celllist',CellList,'celllist',CellList);
     if details.found == 0
         PC.CompareSpikeShape(Clusters, [manualfind(1,:); details.nearest],'tag','CompareShape');
     end
     return;
end
DATA.autolist.listtype = 'autolist'; 

%Step 1. Build initial list by finding good cells then following them
stage = 1;
badid =1;
if isempty(CellList)
    cellq = [];
    [maxq,b] = max(Q(:));
    dcrit = D;
    finddetails.dropcrit = dropcrit;
    dcrit(D < finddetails.dropcrit) = 0;
    dcrit(D >= finddetails.dropcrit) = 1;
    dcrit(DATA.autolist.skipeid,:,:) = 0;
    [maxi, b] = max(I(:).*dcrit(:));
    CellList = zeros(size(Q));
    finddetails.iscrit = iscrit;
    cellid = 1;
    nissue = 0;
    initialdup = [];
    finddetails.fixfit = zeros(size(Q,1),size(Q,2));
    finddetails.fixlist = fixlist;
    [~, nccount] = clust.Check(Clusters,'ncluster')
    finddetails.nc = nccount.nc;


    for j = 1:length(size(fixlist,1))
        e = find(DATA.exptlist == fixlist(j,1));
        if ~isempty(e);
            finddetails.fixfit(e,fixlist(j,2)) = fixlist(j,3);
        end
    end
    while maxi > iscrit(1)
        go = 1;
        [x,y,z] = ind2sub(size(Q),b);
        A = PC.GetClusterInfo(Clusters,[x,y,z]);
%Check that we are on the home probe for this cell
%or that there is a beneift in not being
        [a,b] = max(std(A.MeanSpike.ms'));
        if b ~= y
            B = PC.GetClusterInfo(Clusters,[x,y,0]);
            xc = PC.ShapeCorr(A,B);
            [~,a] = max(xc);
            if CellList(x,b,a) == 0 && I(x,b,a) > 0.8 * I(x,y,z) && I(x,b,a) > iscrit(1) && D(x,b,a) > D(x,y,z)
                y = b;
                z = a;
            else                
                nissue = nissue+1;
                issues(nissue).type = 'isolationprobe';
                issues(nissue).pos = [x,y,z];
                issues(nissue).cellid = cellid;
                                
                fprintf('Starting Cell%d on E%dP%d.%d (I%.2f:D%.2f) But bigger on P%d.%d I%.2f:D%.2f\n',...
                                            cellid,x,y,z,I(x,y,z),D(x,y,z),b,z,I(x,b,a),D(x,b,a));
            end
        end
        
        dupi= [];
        [p,c] = find(squeeze(CellList(x,:,:)) > 0); %current cells on this expt
        for j = 1:length(p)
            B = PC.GetClusterInfo(Clusters,[x p(j) c(j)]);
            if ~iscluster(B)
                fprintf('E%dP%d.%d was a cell, but no cluster!!\n',x,p(j),c(j));
            elseif clust.isdup({A B},'one') %this is already a cell elsewhere
                dupi(j) = I(x,p(j), c(j));
                dupd(j) = D(x,p(j), c(j));
            end
        end
        if ~isempty(dupi) %This is a duplicate of an existing cell, but not yet assigned
            [maxdup, id] = max(dupi);
            badcells(badid).reason = 'CellDup';
            badcells(badid).pos = [x p(id) c(id)];
            if maxdup < I(x,y,z) %why wasnt this square chosen before?  - p(j) c(j) just not a starting point
                nissue = nissue+1;
                issues(nissue).type = 'isolationorder';
                issues(nissue).pos = [x,y,z];
                issues(nissue).cellid = cellid;
            end
            CellList(x,y,z) = -CellList(x,p(id),c(id));
            dcrit(x,p(id),c(id)) = 0;
            initialdup(end+1,:) = [x y z p(id) CellList(x,y,z) cellid];
            go = 0;
        end
        
        if I(x,y,z) > iscrit(1) && go
            fprintf('\n\n%s: Finding Matches for Cell%d on E%dP%dcl%d I%.2fD%.2f\n',datestr(now),cellid,x,y,z,I(x,y,z),D(x,y,z));
            % 'wide' argument means that it also allows for time shifts (caused by triggering).
            % initially this was not used, but seems like a good idea. May make bogus
            % cells made by double triggereing a problem, but should be able to remove
            % those at start.
            oldCellList = CellList;
            [C, Cd,matches, fixdetails] = FindCellSquares(DATA, CellList, Clusters, Q,D,I, [x y z cellid],state,'verbose',finddetails,'corrth',corrth,'wide');
            finddetails.fixfit = fixdetails.fixfit;
            finddetails.fixlist = fixdetails.fixlist;
            if fixdetails.newclusters
                Clusters = getappdata(DATA.toplevel, 'NewClusters');
            end
            CellList(C>0)  = cellid;
            CellList(x,y,z) = cellid;
            a = sub2ind(size(CellList),37, 1, 2);
            cellq(cellid).positions = find(C>0);
            if sum(ismember(a, cellq(cellid).positions))>0 || abs(CellList(a)) > 0
                fprintf('E37P1.2 set here\n');
            elseif abs(CellList(37,1,2)) > 0
                fprintf('E37P1.2 now set\n');
            end
            
            cellq(cellid).Q = Q(x,y,z);
            cellq(cellid).dropi = D(x,y,z);
            cellq(cellid).isolation = I(x,y,z);
            cellq(cellid).probe = y;
            cellq(cellid).expt = x;
            cellq(cellid).cl = z;
            cellq(cellid).matches = matches;
            cellq(cellid).cell = cellid;
            if isempty(Cd) %found matches on all expts
                cellq(cellid).fail(1).id = NaN;
            elseif isfield(Cd{1},'id') && sum(Cd{1}.id) > 0
                cellq(cellid).fail(1).id = Cd{1}.id;
                cellq(cellid).fail(1).xc = max(Cd{1}.pxc);
                cellq(cellid).fail(1).reason = Cd{1}.fail;
            else
                cellq(cellid).fail(1).id = NaN;
            end
            if length(Cd) >1 && isfield(Cd{2},'id') && sum(Cd{2}.id) > 0
                cellq(cellid).fail(2).id = Cd{2}.id;
                cellq(cellid).fail(2).xc = max(Cd{2}.pxc);
                cellq(cellid).fail(2).reason = Cd{2}.fail;
            end
        else
            matches = [];
            go = 0;
        end
        dcrit(find(CellList~=0)) = 0; %This should exlude assigned cells from next lap
        [maxi, b] = max(I(:).*dcrit(:));
        DATA.autolist.CellList = CellList;
        if sum(CellList(:) < 0)
            [maxi, b] = max(I(:).*dcrit(:));
        end
        [~,errs] = celllist.Check(CellList,Clusters);
        if ~isempty(errs)
            id = find(strcmp('cluster',{errs.type}));
            for j = 1:length(id)
                fprintf('!!!!!! CellList and Clusters not matched %s\n',errs(id(j)).s);
            end
            id = find(~strcmp('cluster',{errs.type}));
            if ~isempty(id)
                fprintf('!!!!!! CellList Self Duplicates\n');
            end
        end
        aligndata = {};
        finddetails.cellq = cellq;
        if isfield(errs,'cell')
            id = find([errs.cell] == cellid)
            if ~isempty(id)
                cellq(cellid).listerrs = errs(id);
            end
        end
        LogState(DATA, finddetails);
        if state.watchprogress && go
            
            PC.PlotCellList(DATA,'autolist','lineplot','showfig');
            title(sprintf('%d cells so far',cellid));
            it = findobj(allchild(gcf),'flat','tag','StopButton');
            if isempty(it)
                it = uicontrol(gcf,'style','checkbox','string','stop','units','normalized','position',[0.01 0 0.2 0.1],'Tag','StopButton');
                stoprun = 0;
            else
                stoprun = get(it,'value');
            end
            for j = 1:cellid
                PC.DrawBox(cellq(j).expt,cellq(j).probe,'startpt');
            end
            CellLists{cellid} = CellList;
            drawnow;
            allmatch = cat(2, cellq.matches);
            allpos = cat(1,allmatch.pos);
            na = 1;
            badaligned = 0;
            for j = 1:length(allmatch)
                if isfield(allmatch(j).aligndata,'type') && strcmp(allmatch(j).aligndata.type,'alignfit')
                    aligndata{na} = allmatch(j).aligndata;
                    matchid(na) = j;
                    if aligndata{na}.aligned > 90
                        badaligned(na) = 1;
                    end
                    apos(na,:,:) = aligndata{na}.pos;
                    na = na+1;
                end
            end
            ids = find(allpos(:,1) > 0 & allpos(:,2) > 0);
            ids = sub2ind([size(CellList,1) size(CellList,2)],allpos(ids,1),allpos(ids,2));
            [na,nb] = Counts(ids);
            if max(na) > 1
                fprintf('\n%d locations tested %d times\n',sum(na == max(na)),max(na));
                rpt = nb(na ==max(na));
                [x,y] = ind2sub([size(CellList,1) size(CellList,2)],rpt);
                for j = 1:length(x)
                    fprintf('E%dP%d\n',x(j),y(j));
                end
                fprintf('\n');
            end
            if sum(badaligned)
                fprintf('%d Bad alignments\n',sum(badaligned));
            end
            if stoprun
                fprintf('Stop requested\n');
                yn = gui.Dlg('Stop Requested',DATA.toplevel,{'Stop and Exit' 'Continue'});
            end
        end
        if state.watchprogress > 1 %plot shape consistency matrix
            cellcomps(cellid) = PC.CompareCells(DATA,CellList,Clusters,[cellid cellid],'merged');
            PlotComparison(DATA,cellcomps(cellid),'full');
            drawnow;
        end
        if go
            cellid = cellid+1;
        else
            badcells(badid).matches = matches;
            badcells(cellid).Q = Q(x,y,z);
            badcells(cellid).dropi = D(x,y,z);
            badcells(cellid).isolation = I(x,y,z);
            badcells(cellid).probe = y;
            badcells(cellid).expt = x;
            badcells(cellid).cl = z; 
            badid = badid+1;
        end
        if cellid > stopcell %for testing
            return;
        end
    end
    [~, nccount] = clust.Check(Clusters,'ncluster');
    cid = find(nccount.nc ~= finddetails.nc);
    if ~isempty(cid)
        
    end
    
    id = find(CellList ==0 & ~isnan(I));
    gid = find(dcrit(id) ==0 & D(id) >2 & I(id) > iscrit(2));
    dcrit(D >= 1.8 & CellList == 0) = 1;
    dcrit(DATA.autolist.skipeid,:,:) = 0;
    remI = I(dcrit >0);
    AClusters = getappdata(DATA.toplevel,'AutoClusters')
    finddetails.ChangedClusters{1} = clust.Compare(AClusters, Clusters,'diffonly','noplot');
    finddetails.ReverseClusters{1} = clust.Compare(Clusters, AClusters,'diffonly','noplot');
    setappdata(1,'AutoClusters',Clusters);
    clear AClusters;
    finddetails.missed(1) = sum(remI > 5 & D(dcrit>0) > dropcrit);
    finddetails.missed(2) = sum(remI > 5);
    finddetails.missed(2) = sum(remI > 3);
    [a,b,c] = ind2sub(size(CellList),id(gid));
    
    finddetails.CellLists{1} = CellList;
    finddetails.cellq = cellq;
    finddetails.initial.dup = initialdup;
    buildduplicates =1;
    fixduplicates = 1;
    finddetails.stage = 1;
    finddetails.t(finddetails.stage+1) = now;
    finddetails.initial.aligndata = aligndata;
    finddetails.initial.badcells = badcells;
    setappdata(DATA.toplevel,'AutoCellList',finddetails);
end

if sum(strcmp('save',strargs))
    SaveAutoList(DATA);
end

%Step 2 Find Duplicates. Don't do this piecmeal (expt by expt)
%Look at dupliaces in All Expts, then decide who should replace whom
if finddetails.stage >= stopstage
    return;
end


stage = 2;
if finddetails.stage < 2
    finddetails.celldup = celldup;
    
    if buildduplicates && max(unique(CellList)) > 0
        DATA.selectexpts = [];  %force this to build all expts
        PC.OptionMenu(DATA,CellList,'autoduplicates','nofigures'); %build xcorrCellList and AllDuplicates
    end
    
    if fixduplicates
        [CellList, xdetails] = FixDuplicates(DATA,CellList,D,I,Clusters);
        finddetails.dup.replaced = xdetails.replaced;
        finddetails.dup = CopyFields(finddetails.dup,xdetails,{'replacelist','initialpartiallist' 'warnings'});
%        CellList = FixDuplicates(DATA,CellList,D,I,Clusters);
        CheckDuplicates(DATA, CellList,D,I, Clusters);
    end
    LogState(DATA, finddetails);
    finddetails = SetDetails(DATA, finddetails, CellList, stage, recordstate);
    LogState(DATA, finddetails);
    dup = celllist.FindDuplicates(CellList);
end



%Step 2. Find Duplicates based on
if finddetails.stage >= stopstage
    PC.CheckAlignment(DATA, 'boundary');
    return;
end
if 0 %old stage 2
    %now find duplicates
    [CellList, finddetails.duplicates{1}] = OldFindDuplicates(CellList, I,D);
    finddetails = SetDetails(DATA, finddetails, CellList, 2, recordstate);
end

%Now find squares that look like misclassification/orphans. Start with cells that
%are only one sqares
%
%Can occur in several ways
% E.g 5 matches 4, then 4 ~= 3, or 1. So 5vs never tested.  Then if 3 is a
% cell, it searchs 3+4 = not, then3+5 = match. See lemM301 P6, Expts 3-5

if finddetails.stage >= stopstage
    return;
end

if force.stage == 2 %rerun final check
    CheckDuplicates(DATA, CellList,D,I, Clusters);
end

stage = 3;
if finddetails.stage < 3
    checkcell(ncheck).setcell = 0;
    for cc = 1:2
        [counts, cells] = Counts(CellList(:));
        id = find(counts(:) == cc & cells(:) > 0);
        tmpid = id;
        for j = 1:length(id)
            oldC = CellList;
            setid = find(CellList==cells(id(j)));
            [e,p,cl] = ind2sub(size(CellList),setid);
            %can happen if were two cells each only one expt, then
            %first has two expts on second lap
            if length(e) == 1
                e(2) = e(1);
                p(2) = p(1);
                cl(2) = cl(1);
            else
                [~,a] = min(e);
                [~,b]= max(e);
                e = e([a b]);
                p = p([a b]);
                cl = cl([a b]);
                setid = setid([a b]);
            end
            setcell = 0;
            [x, details{1}] = FindCellForExpt(DATA, e(2)+1, Clusters,Q,D,I,[e(2) p(2) cl(2) cells(id(j))],state,'corrth',0.95,'wide',finddetails,'celllist',CellList);
            if sum(x(:))
                [np,ncl] = find(x > 0);
                setcell = CellList(e(1)+1,np,ncl);
                fprintf('Cell %d (%d,%d,%d) has match on E%dP%d.%d Cell%d\n',cells(id(j)),e(2),p(2),cl(2),e(2)+1,np,ncl,setcell);
                if setcell > 0 && setcell ~= cells(id(j)) %match might not be a cell
                    setcell = -1; %shape match is not cell
                end
            end
            if details{1}.newclusters
                Clusters = getappdata(DATA.toplevel,'NewClusters');
            end
            if ~isfield(details{1},'pxc') %if not in list of goot expts
                details{1}.pxc = [];
            end
            [y, details{2}] = FindCellForExpt(DATA,e(1)-1, Clusters,Q,D,I,[e(1) p(1) cl(1) 0],state,'corrth',0.95,'wide',finddetails,'celllist',CellList);
            if details{2}.newclusters
                Clusters = getappdata(DATA.toplevel,'NewClusters');
            end
            if ~isfield(details{2},'pxc') %if not in list of goot expts
                details{2}.pxc = [];
            end
            if sum(y(:))
                [np,ncl] = find(y > 0);
                if setcell > 0 %if its the same setcell no need to do anything
                    setcell(2) =  CellList(e(1)-1,np,ncl);
                    bothcells = setcell;
                    if setcell(1) ~= setcell(2)
                        if setcell(2) == 0
                            setcell = setcell(1);
                        else
                            nmatch(1) = sum(CellList(:) == setcell(1));
                            nmatch(2) = sum(CellList(:) == setcell(2));
                            fprintf('Cell %d (%d,%d,%d) matches two cells  %d(%d) and %d(%d)!!!',cells(id(j)),e(1),p(1),cl(1),setcell(1),nmatch(1),setcell(2),nmatch(2))
                            if nmatch(1) > nmatch(2) || setcell(2) == 0
                                setcell = setcell(1);
                            else
                                setcell = setcell(2);
                            end
                        end
                        if sum(sum(squeeze(CellList(e,:,:) == setcell))) == 0
                            fprintf(' Using %d\n',setcell);
                            CellList(setid) = setcell;
                            oldcell = unique(abs(CellList(setid)));
                        else
                            xid = find(CellList(e,:,:) == setcell);
                            [a,b,c] = ind2sub(size(CellList(e,:,:)),xid);
                            for k = 1:length(e)
                                if sum(sum(squeeze(CellList(e(k),:,:) == setcell)))
                                    fprintf(' But is on same probe as %d E%d!!!\n',setcell,e(k));
                                else
                                    oldcell = abs(CellList(e(k),p(k),cl(k)));
                                    CellList(e(k),p(k),cl(k)) = setcell;
                                end
                            end
                        end
                    end
                else %setcell(1) was 0. Use 2
                    setcell = CellList(e(1)-1,np,ncl);
                    fprintf('Cell %d (%d,%d,%d) has match on E%dP%d.%d cell%d\n',cells(id(j)),e(1),p(1),cl(1),e(1)-1,np,ncl,setcell);
                    if setcell > 0
                        for k = 1:length(setid)
                            [a,b] = find(squeeze(CellList(e(k),:,:)) == setcell);
                            if isempty(a)
                                CellList(setid(k)) = setcell;
                            else
                                for c = 1:length(a)
                                    fprintf('But Cell%d is already on E%dP%dcl%d',setcell,e(k),a(c),b(c));
                                    PC.CompareSpikeShape(Clusters, [e(k) p(k) cl(k); e(k) a(c) b(c)],'tag','CompareShape');
                                end
                                celldup(end+1).e = e(k);
                                celldup(end).cells = [setcell cells(j)];
                            end
                        end
                    end
                end
                
                if sum(CellList(29,16,:) == 58) > 1
                    fprintf('Cell twice');
                end
                dup = celllist.FindDuplicates(CellList,e,'verbose');
                for k = 1:length(dup)
                    cprintf('red','Cell %d twice at E%d P%dand%d\n',dup(k).cell,dup(k).pos(1:3));
                end
                
            elseif setcell > 0
                CellList(setid) = setcell;
            end
            if setcell == 0
                fprintf('No Matches for Lone Cell %d (%d,%d,%d)\n',cells(id(j)),e(1),p(1),cl(1));
            end
            checkcell(ncheck).cell = cells(id(j));
            checkcell(ncheck).setcell = setcell;
            if size(e,1) > 1
                checkcell(ncheck).elist = e(:);
            end
            checkcell(ncheck).elist = e(:);
            checkcell(ncheck).probelist = p;
            checkcell(ncheck).xc = max([details{1}.pxc details{2}.pxc]);
            checkcell(ncheck).check = 'lone';
            ncheck = ncheck+1;
        end
    end
    finddetails.checkcell = checkcell;


    AClusters = getappdata(DATA.toplevel,'AutoClusters')
    finddetails.ChangedClusters{stage} = clust.Compare(AClusters, Clusters,'diffonly','noplot');
    finddetails.ReverseClusters{stage} = clust.Compare(Clusters, AClusters,'diffonly','noplot');
    setappdata(1,'AutoClusters',Clusters);
    clear AClusters;
    finddetails = SetDetails(DATA, finddetails, CellList, stage, recordstate);
    dup = celllist.FindDuplicates(CellList,[]);
    LogState(DATA, finddetails);
    if sum(strcmp('save',strargs))
        SaveAutoList(DATA);
    end

end


if finddetails.stage < 4 && 0
    %next find places where cell number changes on a probe across a boundary
    %
    [CellList, finddetails.checkcells{4}]  = FixBoundaryChanges(DATA,CellList, Clusters, finddetails);
    dup = celllist.FindDuplicates(CellList);
    finddetails = SetDetails(DATA, finddetails, CellList, 4, recordstate);
    LogState(DATA,finddetails);

end

if finddetails.stage >= stopstage
    return;
end
if finddetails.stage < 5
    %next find places where cell number changes on a probe across a boundary
    %CheckCellShapes builds correlation matrices. ShapeMerge does the
    %merges identified in finddetails.comp
    if 1
        [CellList, finddetails]  = CheckCellShapes(DATA, CellList, Clusters, finddetails,'noplot');
    end
    CellList = ShapeMerge(DATA, CellList, Clusters, finddetails);
    dup = celllist.FindDuplicates(CellList);
    if isappdata(DATA.toplevel,'SelfDetails')
        x = getappdata(DATA.toplevel,'SelfDetails');
    else
        x = CheckSelfShapes(DATA, CellList,Clusters, finddetails);
        setappdata(DATA.toplevel,'SelfDetails',x);
    end
    if state.watchprogress > 1
        PlotComparison(DATA,x([x.merge] > 1),'fullwait');
    end
    finddetails = SetDetails(DATA, finddetails, CellList, 5, recordstate);
    LogState(DATA, finddetails);
    if sum(strcmp('save',strargs))
        SaveAutoList(DATA);
    end
end

if finddetails.stage >= stopstage
    return;
end
if finddetails.stage < 6
    %next find places where cell number changes on a probe across a boundary
    %
    %[CellList, finddetails]  = CompareCellShape(CellList, Clusters, finddetails);
   % finddetails = SetDetails(DATA, finddetails, CellList, 6, recordstate);
end

if finddetails.stage >= stopstage
    return;
end
if finddetails.stage < 7
    %Now remove cells that appear too briefly
    [a,b] = Counts(CellList);
    a = a(b>0);
    b = b(b>0);
    bid = b(a < minexpts);
    clearcell = [];
    nx = 1;
    for j = 1:length(bid)
        id = find(CellList == bid(j));
        [e,p,cl] = ind2sub(size(CellList),id);
        fprintf('Deleting cell %d (E%dP%d): Only %d Expts\n',bid(j),e(1),p(1),length(e));
        clearcell(nx).cell = bid(j);
        clearcell(nx).cleared = id;
        clearcell(nx).eid = e;
        clearcell(nx).pid = p;
        CellList(id) = 0;
        nx = nx+1;
    end
    
    finddetails.clearcell = clearcell;
    finddetails = SetDetails(DATA, finddetails, CellList, 7,recordstate);
    LogState(DATA, finddetails);
    DATA.autolist.CellList = CellList;
    dup = celllist.FindDuplicates(CellList);
    LogState(DATA, finddetails);
end



if finddetails.stage >= stopstage
    return;
end
if finddetails.stage < 9
        [CellList, finddetails.fixboundary] = FixBoundaryChanges(DATA,CellList, Clusters, finddetails);
        finddetails = SetDetails(DATA, finddetails, CellList, 9, recordstate);
        LogState(DATA, finddetails);
end

if finddetails.stage >= stopstage
    return;
end
finddetails.state = state;
if finddetails.stage < 10    
        [CellList, finddetails.fixgaps] = FixGaps(DATA,CellList, Clusters, Q,D,I, finddetails);
        finddetails = SetDetails(DATA, finddetails, CellList, 10, recordstate);
        LogState(DATA, finddetails);
end

%finddetails.StartClusters = Clusters;

if finddetails.stage >= stopstage
    return;
end

if force.swapclusters >= 0
    if finddetails.stage < 11
        [CellList, finddetails.cellid] = ReNumberCellList(CellList);
        [CellList, details] = celllist.Renumber(CellList,'nloop',10000); %keep close #s far apart
        finddetails.cellid = finddetails.cellid(details.cellid);
        if alignclusters %only do this when everything else workds
            [CellList, Clusters,finddetails.clusterswaps] = AlignClustersForCell(CellList, Clusters,D);
        end
        %re-calc params after swapping clusters around
        [Q, D, I] = PC.MakeVarMatrix(DATA,Clusters,'cellquality','dropi','defaultisolation');
        DATA.autolist.dropi = D;
        
        %run again to check
        AlignClustersForCell(CellList, Clusters,D,'noswap');
        
        [CellList, finddetails.finaldup] = FixDuplicateSquares(DATA,CellList, Clusters,finddetails);
        finddetails.Clusters = Clusters;
        finddetails = SetDetails(DATA, finddetails, CellList, 11, recordstate);
        LogState(DATA, finddetails);        
    end
end


if force.addmu >= 0 && finddetails.stage < 11
    DATA.autolist.muCellList = FindMuCells(DATA, CellList, Clusters, Q, D,I, state);
    finddetails.muCellList = DATA.autolist.muCellList;
    setappdata(finddetails.toplevel,'AutoCellList',finddetails);
end

finddetails.CellList = CellList;
save([DATA.datadir '/PreCheckClusters.mat'],'Clusters');
if sum(strcmp('save',strargs))
    SaveAutoList(DATA);
end

if finddetails.stage >= stopstage
    return;
end
state.watchprogress = 1;
finddetails.stage = 11;
finddetails = PC.CheckAlignment(DATA, finddetails,'apply');

stage = 11;
if isappdata(DATA.toplevel,'NewClusters');
    Clusters = getappdata(DATA.toplevel,'NewClusters');    
    AClusters = getappdata(DATA.toplevel,'AutoClusters')
    finddetails.ChangedClusters{fstage} = clust.Compare(AClusters, Clusters,'diffonly','noplot');
    finddetails.ReverseClusters{stage} = clust.Compare(Clusters, AClusters,'diffonly','noplot');
    setappdata(DATA.toplevel,'AutoClusters',Clusters);
    clear AClusters;
end
finddetails = SetDetails(DATA, finddetails, CellList, 11, recordstate);
LogState(DATA, finddetails);
finddetails = AddChecks(DATA, finddetails);

%check again for possibel duplicates
%if stage = 14 uses save duplicatee list from previous attemp
%stage 12 recalculates
if ismember(force.stage,[13 14])%check again for duplicates
    if force.stage == 14
        DATA.selectexpts = [];  %force this to build all expts
        PC.OptionMenu(DATA,CellList,'autoduplicates','nofigures'); %build xcorrCellList and AllDuplicates
    end
    [CellList, xdetails] = FixDuplicates(DATA,CellList,D,I,Clusters);
    finddetails.bdup.replaced = xdetails.replaced;
end
DATA.autolist.CellList = CellList;
finddetails = SetDetails(DATA, finddetails, CellList, 12, recordstate);
DATA.autolist = CopyFields(DATA.autolist,finddetails,{'cellid' 'InitialList' 'checkcell' 'clearcell'});

PC.CheckAutoList(DATA);

if sum(strcmp('save',strargs))
    SaveAutoList(DATA);
end
if nargout == 0
    SetData(DATA);
end



function FX = AddChecks(DATA, FX)
 
 CellList = FX.CellList;
 D = FX.D;
 I = FX.I;
cellA = FX.A; 
    
FX.candidates = find(D > 1.7 & I > 2.5 & CellList ==0); 
[e,p,cl] = ind2sub(size(CellList),FX.candidates);
id = find(ismember(e,DATA.autolist.eid));
FX.candidates = FX.candidates(id);


cid = find(D < 1.5 & I > 3.5 & CellList ==0 & cellA > 0.9); 
[e,p,cl] = ind2sub(size(CellList),cid);
id = find(ismember(e,DATA.autolist.eid));
FX.dropcells = cid(id);

FX.selfdup = [];
for e = 1:size(CellList,1)
    [nc, cells] = Counts(abs(CellList(e,:,:)));
    id = find(nc > 1 & cells > 0);
    for k = 1:length(id)
        [ei,p,cl] = celllist.find(abs(CellList),cells(id(k)),'expts',e);
        [a,b] = Counts(p);
        if max(a) > 1
            [cp, ccl] = find(squeeze(CellList(e,:,:)) == cells(id(k))); %true cell
            aid = find(a > 1);
            for c = 1:length(aid)
                pcl = find(squeeze(CellList(e,b(aid(c),:)) == -cells(id(k)))); %on this probe
                [dp, dcl] = find(squeeze(abs(CellList(e,:,:))) == cells(id(k))); %dups
                [dp, ix] = setdiff(p,b(aid(c))); %dup  or cellon another probe
                if isempty(dp)
                    [dp, dcl] = find(squeeze(CellList(e,:,:)) == -cells(id(k))); %dups
                    FX.selfdup(end+1,:) = [e b(aid(c)) dcl(1) cells(id(k)) cp ccl];
                else
                    FX.selfdup(end+1,:) = [e b(aid(c)) dcl(1) cells(id(k)) dp(1) dcl(ix(1))];
                end
            end
        end
    end
end

FX.ratecheck = PC.CheckAllRateSequences(DATA);

    function [CellList, fixes] = FixGaps(DATA, CellListIn, Clusters, Q, D, I, finddetails)
%Once cellist has been filled, find gaps. can afford a lower threshodl
%now for xcorr, as will not steal a cell, only take empty cells

      fixes(1).cell = 0;
      fixes(1).pos = 0;
      CellList = CellListIn;
      [counts, cells] = Counts(CellList(:));
      cid = find(counts(:) > 3 & cells(:) > 0);
      for j = 1:length(cid)
          oldC = CellList;
          cell = cells(cid(j));
          [expts,p,cl] = celllist.find(CellList,cell);
          de = find(diff(expts) > 1);
          dp = p(de);
          dcl = cl(de);
          de = expts(de);
          xde = de+1;
          for e = 1:length(de)
              [CellList, fixes(end+1)] = FixGap(DATA,CellList,Clusters, Q,D,I, xde(e), [de(e) dp(e) dcl(e) cell],finddetails);
          end
          e = expts(1);
          while e > 1
              [CellList, fixes(end+1)] = FixGap(DATA,CellList,Clusters, Q,D,I, e-1, [e p(1) cl(1) cell],finddetails);
              if fixes(end).cell == 0
                  e = 0;
              else
                  e = e-1;
              end
          end
          while e < size(CellList,1)-1
              [CellList, fixes(end+1)] = FixGap(DATA,CellList,Clusters, Q,D,I, e, [e+1 p(1) cl(1) cell], finddetails);
              if fixes(end).cell == 0
                  e = NaN;
              else
                  e = e+1;
              end
          end
      end
      fixes = fixes([fixes.cell] > 0);

function SaveAutoList(DATA)
    finddetails = RemoveHandles(getappdata(DATA.toplevel,'AutoCellList'));
    savename = sprintf('%s/AutoCellList.mat',DATA.autodatadir);
    if isfield(finddetails,'Clusters')
        Clusters = finddetails.Clusters;
    elseif isappdata(DATA.toplevel,'NewClusters')
        Clusters = getappdata(DATA.toplevel,'NewClusters')
    else
        Clusters = PC.GetValue(DATA, 'Clusters');
    end
    if ~isfield(finddetails,'clusterswaps')
        finddetails.clusterswaps = ListClusterSwaps(Clusters);
    end
    finddetails = rmfields(finddetails,{'Clusters','StartClusters'});
    finddetails.savetime = now;
    fprintf('Saving %s\n',savename);
    save(savename,'-v7.3','finddetails');
    
 function finddetails = SetDetails(DATA, finddetails, celllist, stage, rec)
 %finddetails = SetDetails(DATA, finddetails, celllist, stage, rec)
    finddetails.CellLists{stage} = celllist;
    DATA.autolist.CellList = celllist;
    finddetails.stage = stage;
    finddetails.CellList = celllist;
    finddetails.t(finddetails.stage+1) = now;
    if rec >= stage
        setappdata(finddetails.toplevel,'AutoCellList',finddetails);
    end

    
    
function DATA = FixList(DATA, Clusters, CellList, stage, varargin)
    
  finddetails = getappdata(DATA.toplevel,'AutoCellList');
  if stage == 100
%step 1 just builds a list of expts/probes where there are duplicates (including any already detected)      
  [a,b] = celllist.FixDuplicates(CellList,Clusters, [],'all','verbose');
%step2 determines appropriate fixes
  [~, fixes] = celllist.FixDuplicates(CellList,Clusters, b);  
%make backup of Clusters, and then apply fixes
  if ~isappdata(DATA.toplevel,'Clusterbackup')
      setappdata(DATA.toplevel,'Clusterbackup',Clusters);
  end
  [C, Cl] = celllist.FixDuplicates(CellList,Clusters, fixes);
  
  setappdata(DATA.toplevel,'AutoClusters',Cl);
  finddetails.CellList = C;
  finddetails.CellLists{12} = C;
  DATA.autolist.CellList = C;
  end
  
  setappdata(DATA.toplevel,'AutoCellList',finddetails);
  PC.PlotCellList(DATA,'autolist','showfig');
  if nargout == 0
      SetData(DATA);
  end

  
  
  
function [CellList, fix] = FixGap(DATA, CellList, Clusters, Q, D, I, e, tid, finddetails)
    fix.cell = 0;
    fix.pos = 0;
    celln = tid(end);
    
    [x, details] = FindCellForExpt(DATA,e, Clusters,Q,D,I,tid,finddetails.state,'corrth',0.95,'wide',finddetails,'celllist',CellList);
    if details.found
        c = CellList(details.id(1),details.id(2),details.id(3));
        if c == 0
            fix.cell = celln;
            fix.pos = details.id;
            CellList(details.id(1),details.id(2),details.id(3)) = celln;
        else
            fprintf('Cell %d fills gap in cell%d at E%dP%dcl%d: xc %.3f\n',c,celln,e,tid(2:3),max(details.pxc));
        end        

    end


function [CellList, checkcell] = FixBoundaryChanges(DATA, CellList, Clusters, finddetails)
    corrth = finddetails.corrth;
    ncheck = 1

    X.toplevel = gcf;
    nc = 0;
    for p = 1:size(CellList,2)
        oldC = CellList;
        C = squeeze(CellList(:,p,:));
        [counts, cells] = Counts(C);
        counts = counts(cells > 0);
        cells = cells(cells > 0);
        replaced = zeros(1,max(cells));
        if sum(counts > 1) > 1
            for k = 2:length(cells)
                for j = 1:k-1
                    cellj = cells(j);
                    if replaced(cellj)
                        cellj = replaced(cellj);
                    end
                    cellk = cells(k);
                    if replaced(cellk)
                        cellk = replaced(cellk);
                    end
                    [pa,ca] = find(C == cellj);
                    [pb,cb] = find(C == cellk);
                    dp = bsxfun(@minus,pa,pb');
                    [a,b] = find(abs(dp)<2);
                    nsame = sum(dp(:) ==0);
                    if length(a) > nsame * 2
                        matched = 0;
                        for c = 1:length(a) %check each transiton
                            if sum(C(pa(a(c)),:) == cellk) == 0 && sum(C(pb(b(c)),:) == cellj) == 0 %if both cells are on the same probe
                                fprintf('P%d Cell label changes E%d cell%d(cl%d)->E%d cell%d(cl%d)\n',p,pa(a(c)),cellj,ca(a(c)),pb(b(c)),cellk,cb(b(c)));
                                % only check xcorr for this cell
                                A = PC.GetClusterInfo(Clusters,pb(b(c)),p,cb(b(c)),'allexpt');
                                B = PC.GetClusterInfo(Clusters,pa(a(c)),p,ca(a(c)),'allexpt');
                                %allow for shift in triggering at this point
                                [xc, details{1}] = PC.ShapeCorr(A,B,'delays',5,'proberange',finddetails.corrproberange,'drifts',1,finddetails.ArrayConfig);
                                checkcell(ncheck).cell = cellk;
                                checkcell(ncheck).setcell = cellj;
                                checkcell(ncheck).elist = [pb(b(c)) pa(a(c))];
                                checkcell(ncheck).probelist = [p p];
                                checkcell(ncheck).xc = xc;
                                checkcell(ncheck).details = details{1};
                                checkcell(ncheck).check = 'boundary';
                                nc = nc+1;
                                cellcomps(nc) = CompareCells(DATA,CellList,Clusters,[cellj cellk]);
                                cellcomps(nc).p = p;

                                if matched
                                    fprintf('Already Set Cell %d -> %d\n',cellk,cellj);
                                elseif cellcomps(nc).merge > 0
                                    fprintf('Setting Cell %d -> %d\n',cellk,cellj);
                                    replaced(cellk) = cellj;
                                    xid = sub2ind(size(CellList),pb, ones(size(pb)) .* p, cb);
                                    for e = 1:length(pb)
                                        if sum(sum(CellList(pb(e),:,:) == cellj)) == 0
                                            CellList(pb(e),p,cb(e)) = cellj;
                                            did = find(CellList(pb(e),:,:) == -cellk);
                                            CellList(did) = -cellj;
                                            matched = 1;
                                        else
                                            fprintf('Cell%d already defined on E%dP%d. Clusters %d and %d same shape?\n',cellj,pb(b(c)),p,ca(a(c)),cb(b(c)));
                                        end
                                        if sum(sum(CellList(pb(e),:,:)) == cellj) > 1
                                            fprintf('Duplicated celld %d\n',cellj);
                                        end
                                    end
%make C reflect new assignements                                    
                                    C = squeeze(CellList(:,p,:));
                                    checkcell(ncheck).replacelist = pb;
                                    xid = find(CellList(:) == cellk);
                                    c = length(a); %break out of loop now cell has been switched
                                    if ~isempty(xid)
                                        fprintf('Cell %d is still present elsewoere\n',cellk);
                                    end
                                    checkcell(ncheck).unreplacelist = xid;
                                elseif xc > corrth(2)
                                    fprintf('Cells %d -> %d match at boundary %d/%d, but not elsewhere\n',cellk,cellj,pa(a(c)),pb(b(c)));
                                elseif isnan(xc) %missing some data
                                    if isempty(B)
                                        fprintf('No Data for E%dP%dcl%d\n',pa(a(c)),p,ca(a(c)));
                                    end
                                else
                                    fprintf('Cell %dE%dP%d and Cell %dE%dP%d seem different (xc%.3f)\n',cellj,pa(a(c)),p,cellk,pb(b(c)),p,xc);
                                    GetFigure('CompareShape');
                                    CompareSpikeShape(A,B);
                                    checkcell(ncheck).replacelist = 0;
                                end
                                ncheck = ncheck+1;
                            elseif finddetails.verbose
                                fprintf('E%d Cells %d,%d on Probe %d\n',pa(a(c)),cellj,cellk,p);
                            end
                        end
                        np = 0;
                        e = celllist.find(CellList,cells(j));
                        if max(Counts(e)) > 1
                            fprintf('Duplicated cell %d\n',cellj);
                        end
                    end
                end
            end
        end
        C = squeeze(CellList(:,p,:));
        [counts, cells] = Counts(C);
        counts = counts(cells > 0);
        cells = cells(cells > 0);
        for e = 1:size(C,1)
            [counts, cells] = Counts(C(e,:));
            counts = counts(cells > 0);
            if max(counts) > 1
                cid = find(counts> 1);
                fprintf('Cell Duplicated on E%dP%d: %d\n',e,p,cells(cid));
            end
        end
    end

    checkcell;
 
 function [CellList, finddetails] = ShapeMerge(DATA, CellList, Clusters, finddetails)
     comp = finddetails.cellcomps;
     id = find([comp.merge] > 0);
     comp = comp(id);
     mergegroup{1} = comp(1).cells;
     for j = 1:length(comp)
         found = 0;
         for k = 1:length(mergegroup)
             if sum(ismember(comp(j).cells,mergegroup{k}))
                 found = found+1;
                 mergegroup{k} = union(mergegroup{k},comp(j).cells);
             end
         end
         if ~found
             mergegroup{end+1} = comp(j).cells;
         end
     end
     for j = 1:length(mergegroup)
         [e,p,cl] = celllist.find(CellList,mergegroup{j});
         [a,b] = Counts(e);
         if max(a) > 1
         else
             fprintf('Merging Cells %s\n',sprintf(' %d',mergegroup{j}));
             id = sub2ind(size(CellList),e,p,cl);
             CellList(id) = mergegroup{j}(1);
             [e,p,cl] = celllist.find(CellList,mergegroup{j}(1));
             [a,b] = Counts(e);
             if max(a) > 1
                 fprintf('Duplicated Cell\n');
             end
         end
     end
     
function [CellList, finddetails] = CheckCellShapes(DATA, CellList,  Clusters, finddetails,varargin)

    plottype = 'default';
    j = 1;
    while j <= length(varargin)
       if strcmpi(varargin{j},'noplot')
            plottype = 'none';
        end
        j = j+1;
    end
    
    nc = 1;
    for p = 1:size(CellList,2)
        C = squeeze(CellList(:,p,:));
        [counts, cells] = Counts(C);
        counts = counts(cells > 0);
        cells = cells(cells > 0);
        if sum(counts > 1) > 1
            for k = 2:length(cells)
                for j = 1:k-1
                    cellcomps(nc) = CompareCells(DATA,CellList,Clusters,cells([j k]),varargin{:});
                    cellcomps(nc).p = p;
                    nc = nc+1;
                end
            end
        end
    end
    finddetails.cellcomps = cellcomps;
    
function cellcomps = CheckSelfShapes(DATA, CellList,  Clusters, finddetails)
    nc = 1;
    cells = unique(CellList);
    cells = cells(cells>0);
    for j = 1:length(cells)
        if sum(CellList(:) == cells(j)) > 1
            fprintf('Checking shape consistency for cell %d\n',cells(j));
            cellcomps(nc) = CompareCells(DATA,CellList,Clusters,cells([j j]),'merged');
            nc = nc+1;
        end
    end
    finddetails.cellcomps = cellcomps;
    
    
function mykey(a,b);
     type = get(gcf,'SelectionType');
     if strcmp(b.Key,'shift')
         fprintf('Jumping to the End\n')
         F = GetFigure(a);
         setappdata(F,'waitforbutton',0);
         uiresume(F);
     end
     
     
function AddNextButton(F)    

    
    it = findobj(allchild(F),'flat','tag','NextButton');
    if isempty(it)
        it = uicontrol(F,'units','normalized','position',[0 0.9, 0.2 0.1],'tag','NextButton','string','next','callback','uiresume(gcbf)');
        set(it,'keypressfcn',@mykey);
        setappdata(F,'waitforbutton',1);
    end

function F = PlotComparison(DATA, comp,type)
    if strcmp(type,'none')
        F = 0;
        return;
    end
    if length(comp) > 1
        for j = 1:length(comp)
            if ~strcmp(type,'none')
            F = PlotComparison(DATA,comp(j),type);
            end
            if strcmp(type,'fullwait')
                if j == 1
                    setappdata(F,'waitforbutton',1);
                end
                AddNextButton(F);
                a = getappdata(F,'waitforbutton');
                if a
                    uiwait(gcf);
                else
                    fprintf('Got Jump\n');
                    type = 'none';
                    break;
                end
            end
        end
        return;
    end
    F = GetFigure('CompareShapeMatrix','parent',DATA.toplevel);
    hold off; 
    h = imagesc(comp.xc);
    minstep = min(diag(comp.xc,1));
    stepstr = sprintf('Worst Step %.3f',minstep);
    set(gca,'ytick',1:length(comp.ea),'yticklabel',num2str(comp.ea));
    set(gca,'xtick',1:length(comp.eb),'xticklabel',num2str(comp.eb));
    caxis([max([min(comp.xc(:)) 0.5]) 1]);
    colorbar;
    xlabel(sprintf('Cell%dP%d',comp.cells(1),prctile(comp.pa,50)));
    ylabel(sprintf('Cell%dP%d',comp.cells(2),prctile(comp.pb,50)));
    if strncmp(type,'full',4)
        setappdata(gcf, 'CompData',comp);
        set(h,'buttondownfcn',{@HitCompIm});
    end
    if comp.merge > 0
        if isempty(comp.nonmatch)
            title(sprintf('Merge %s',stepstr));
        else
            title(['Merge but check E' sprintf('%d ',comp.nonmatch) stepstr]);
        end
    else
        if isempty(comp.goodmatch)
            title(sprintf('NoMerge %s',stepstr));        
        else
            title(['NoMerge But Check'  sprintf('%d ',comp.goodmatch) stepstr]);        
        end            
    end
        
    
function HitCompIm(src,b)

    DATA = GetDataFromFig(src);
    Clusters = getappdata(DATA.toplevel,'AutoClusters');
    finddetails = getalldata(DATA.toplevel,'AutoCellList');
    xy = get(gca,'currentpoint');
    comp = getappdata(gcf,'CompData');
    a = round(xy(1,2));
    b = round(xy(1,1));
    A = PC.GetClusterInfo(Clusters,[comp.ea(a) comp.pa(a) comp.ca(a)]);
    B = PC.GetClusterInfo(Clusters,[comp.eb(b) comp.pb(b) comp.cb(b)]);
    if comp.ea(a) == comp.eb(b)
        e = CalcEfficacy(A.times,B.times);
        xc = xcorrtimes(A.times,B.times);
        GetFigure('xcorr','parent',DATA.toplevel);
        plot(xc);
        title(sprintf('E%d',comp.ea(a)));
    else
        GetFigure('CompareShape','parent',DATA.toplevel);
        PC.CompareSpikeShape({A B});
        [xc, details] = PC.ShapeCorr(A,B, DATA.ArrayConfig,'shifts',5,'drifts',1,'proberange',finddetails.corrproberange);
        title(sprintf('Best xc %.3f',xc));
    end
    

    
    
function comp = CompareCells(DATA, CellList, Clusters, cid, varargin) 
    
    plottype = '';
    checkmerge = 0;
    if isfield(DATA,'ArrayConfig')
        Array = DATA.ArrayConfig;
    else 
        Array = [];
    end
    j = 1;
    while j <= length(varargin)
        if strncmp(varargin{j},'merged',4);
           checkmerge = 1;
        elseif strcmpi(varargin{j},'noplot')
            plottype = 'none';
        elseif strncmp(varargin{j},'plot',4);
           plottype = 'full';
        end
        j = j+1;
    end
    [ea, pa, ca] = celllist.find(CellList,cid(1));
    [eb, pb, cb] = celllist.find(CellList,cid(2));
    xc = [];
    efficacy = [];
    for j = 1:length(ea)
        for k = 1:length(eb)
            C{1} = PC.GetClusterInfo(Clusters,[ea(j) pa(j) ca(j)],'allexpt');
            C{2} = PC.GetClusterInfo(Clusters,[eb(k) pb(k) cb(k)],'allexpt');
            [xc(j,k), details] = PC.ShapeCorr(C{1},C{2},'delays',5,'drifts',2,Array,'proberange',DATA.autolist.corrproberange);
            if xc(j,k) < 0.8
            end
            if ea(j) == eb(k)  && checkmerge == 0%same expt, check efficacy
                e = CalcEfficacy(C{1},C{2});
                efficacy(ea(j)) = max(e(1:2));
            end
        end
    end
    comp.p = prctile(cat(1,pa, pb),50);
    comp.cells = cid;
    comp.ea = ea;
    comp.eb = eb;
    comp.pa = pa;
    comp.pb =pb;
    comp.ca = ca;
    comp.cb = cb;
    comp.xc = xc;
    comp.nonmatch = [];
    comp.goodmatch = [];
    comp.merge = -1;
    comp.efficacy = efficacy;
    eboth = intersect(ea,eb);
    gid = find(xc > 0.95);
    bid = find(xc < 0.8);
    xcm = mean(xc(:));
    if checkmerge
        comp.merge = 1;
        xca = mean(xc,2);
        id = find(xca < 0.9);
        comp.nonmatch = ea(id);
        if ~isempty(id)
            comp.merge = 2;
        end
        if size(xc,1) > 1
            transitions = diag(xc,1);
            id = find(transitions < 0.9);
            comp.badstep = ea(id);
        else
            comp.badstep = [];
        end
%if a whole expt it bad, guarantees that some transitions will be        
        if ~isempty(id) && comp.merge == 0
            comp.merge = 3;
        end
        return;
    end
    if ~isempty(gid)
        did = find(efficacy(eboth) < 0.1);
        sid = find(efficacy(eboth) > 0.2);
        for j = 1:length(eboth)
            
        end
        if ~isempty(did)
            if length(sid) > length(did)
                fprintf('Cells %d and %d xc%.3f could merge %.3f %d/%d sync, %d not\n',cid,xcm,length(sid),length(eboth),length(did));
                comp.merge = 3
            else
                fprintf('Cells %d and %d: xc%.3f. %d/%d efficacy < 0.1, %d > 0.2\n',cid,xcm,length(did),length(eboth),length(sid));
                comp.merge = -1;
                if ~isempty(sid) %need to check wehre it appeas to be the same; 
                    comp.merge = -2;
                end
            end
        elseif ~isempty(bid)
            fprintf('Cells %d and %d could merge %.3f %d good, but %d bad',cid,xcm,length(gid),length(bid));
            [a,b] = ind2sub(size(xc),bid);
            nx = [length(unique(a)) length(unique(b))];
            if nx(1) == 1 && size(xc,1) > 2
                fprintf('E%dP%d cell %d is problem',ea(unique(a)),pa(unique(a)),cid(1));
            end
            if nx(2) == 1 && size(xc,2) > 2
                fprintf('E%dP%d cell %d is problem',eb(unique(b)),pb(unique(b)),cid(2));
            end
            fprintf('\n');
            xca = mean(xc,2);
            id = find(xca < 0.9);
            if ~isempty(id)
                fprintf('Cell%d bad at %s',cid(1),sprintf('%d ',ea(id)));
                comp.nonmatch = [comp.nonmatch ea(id)'];
            end
            id = find(xca > 0.95);
            if ~isempty(id)
                comp.goodmatch = [comp.goodmatch ea(id)'];
            end
            xcb = mean(xc,1);
            id = find(xcb < 0.9);
            if ~isempty(id)
                fprintf('Cell%d bad at %s',cid(2),sprintf('%d ',eb(id)));
                comp.nonmatch = [comp.nonmatch eb(id)'];
            end
            id = find(xcb > 0.95);
            if ~isempty(id)
                comp.goodmatch = [comp.goodmatch eb(id)'];
            end
            if xcm > 0.9 || prctile(xc(:),20) > 0.9
                fprintf('Should Merge (%s)\n',datestr(now));
                comp.merge = 2;
                if ~isempty(comp.nonmatch)
                    comp.merge = 4;
                end
            else
                comp.merge = 0;
            end
        else
            fprintf('Cells %d and%d should merge mean xc: %.2f\n',cid,xcm);
            comp.merge = 1;
        end
    else
        comp.merge = 0;
    end
    if comp.merge <= 0 && ~isempty(comp.goodmatch)
        comp.merge = -2;
    end
    PlotComparison(DATA,comp,plottype);
    drawnow; 
    
function [C, Clusters, clusterswaps] = AlignClustersForCell(C, Clusters, D, varargin)

    noswap = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'noswap',6)
            noswap = 1;
        end            
        j = j+1;
    end
    
    cells = unique(C(C>0));
    swaps = [];
    runtwice = 0;
    for j = 1:length(cells)
        cid = find(C==cells(j));
        [a,b,c] = ind2sub(size(C),cid);
        cl = mode(c);
        p = mode(b);
        if cl > 2
            nset(1) = sum(C(a,p,1)>0 & C(a,p,1)~= cells(j));
            if nset(1) == 0
                cl = 1;
            else
                nset(2) = sum(C(a,p,2)>0  & C(a,p,2)~= cells(j));
                if nset(1) == 0
                    cl = 2;
                end
            end
        end
        if cl ~= mode(c)
            fprintf('Using cl %d for cell %d\n',cl,cells(j))
        end
        for k = 1:length(cid)
%could restrict swapts to cases where this probe == p?              
            if c(k) ~= cl && b(k) == p
                if noswap
                    fprintf('Cell%d Want to Swap E%dP%d %d <-> %d',cells(j),a(k),b(k),cl,c(k));
                    ok = 0;
                else
                    fprintf('Cell%d Swapping E%dP%d %d <-> %d',cells(j),a(k),b(k),cl,c(k));
                    [Clusters{a(k)}{b(k)}, ok] = clust.Swap(Clusters{a(k)}{b(k)},cl,c(k));
                end
                if ok
                    swaps(end+1,:) = [a(k) b(k) cl c(k)];
                    C(a(k),b(k),c(k)) = C(a(k),b(k),cl);
                    C(a(k),b(k),cl) = cells(j);
                    d = D(a(k),b(k),c(k));
                    D(a(k),b(k),c(k)) = D(a(k),b(k),cl);
                    D(a(k),b(k),cl) = d;
                    fprintf('...OK\n');
                else
                    fprintf('...Cancelled\n');
                end
            end
        end
    end

if nargout < 3
    return;
end
if isempty(swaps)
    clusterswaps = {};
else
    clusterswaps{max(swaps(:,1))}{max(swaps(:,2))} = [];
    for j = 1:size(swaps,1)
        if length(clusterswaps{swaps(j,1)}) < swaps(j,2) || isempty(clusterswaps{swaps(j,1)}{swaps(j,2)})
            clusterswaps{swaps(j,1)}{swaps(j,2)} = swaps(j,3:end);
        else
            clusterswaps{swaps(j,1)}{swaps(j,2)}(end+1,:) = swaps(j,3:end);
        end
    end
end
%dont see why the second run through is there.  Causes double swaps which
%are a problem.  So leave off until see need.
if runtwice    
    for j = 1:length(cells)
        cid = find(C==cells(j));
        [a,b,c] = ind2sub(size(C),cid);
        cl = mode(c);
        for k = 1:length(cid)
            if c(k) ~= cl
                fprintf('Cell%d Swapping E%dP%d %d <-> %d\n',cells(j),a(k),b(k),cl,c(k));
                [Clusters{a(k)}{b(k)}, ok] = clust.Swap(Clusters{a(k)}{b(k)},cl,c(k));
                if ok
                    swaps(end+1,:) = [a(k) b(k) cl c(k)];
                    C(a(k),b(k),c(k)) = C(a(k),b(k),cl);
                    C(a(k),b(k),cl) = cells(j);
                    d = D(a(k),b(k),c(k));
                    D(a(k),b(k),c(k)) = D(a(k),b(k),cl);
                    D(a(k),b(k),cl) = d;
                else
                    
                end
            end
        end
    end
end

function [C, Cd, matches, fixdetails] = FindCellSquares(DATA,CellList, Clusters, Q,D,I,id,state, varargin)

%Clusters = getappdata(DATA.toplevel,'Clusters');
C = zeros(size(Q));
C(id(1),id(2),id(3)) = 1;
Cd ={};
matches = [];
cellid = id(4);
verbose = state.verbose;

j = 1;
while j <= length(varargin)
    if isfield(varargin{j}, 'ArrayConfig')
        FX = varargin{j};
        ArrayConfig = varargin{j}.ArrayConfig;
        fixfit = FX.fixfit;
    elseif isfield(varargin{j},'xx')
        FX = varargin{j};
    end
    j = j+1;
end
if isfield(FX,'fixlist')
    fixdetails.fixlist = FX.fixlist;
else
    fixdetails.fixlist = [];
end
fixdetails.fixfit = fixfit;
fixdetails.newclusters = 0;
tid = id;
goodid = tid;
goode = id(1);
bade = 0;
ne = 0;
nbad = 0;
tC = CellList;
for e = id(1)+1:size(Q,1);
    ne = ne+1;
    cid = find(C>0);
    tC(cid) = C(cid);
    if ismember(e,DATA.autolist.eid)
        chkdetails = [];
        [x, details] = FindCellForExpt(DATA, e, Clusters,Q,D,I,goodid,state,varargin{:},DATA,'celllist',CellList);
        if details.newclusters
            Clusters = getappdata(DATA.toplevel,'NewClusters');
            fixdetails.newclusters = fixdetails.newclusters+1;
        end
        if details.fixfit
            fixdetails.fixlist(end+1,:) = [details.id(1:2) details.fixfit 2 now];
            fixdetails.fixfit(details.id(1),details.id(2)) = details.fixfit;            
        end
        %do NOT overwrite existing cells wiht new ones - process starts
        %with the best
        oldx = squeeze(CellList(e,:,:));
        [ap,acl] = find(x >0);
        x(abs(oldx)>0) = 0;
        C(e,:,:) = x;
        tC(e,:,:) = x .* cellid;
        if details.found
            oldid = goodid;
            [p,cl] = find(x >0);
            if isempty(p)
                p = ap(1);
                cl = acl(1);
                matches(ne).good = 2; %duplicate
            else                
                matches(ne).good = 1;
                goodid = [e p cl cellid];
                goode = e;
                [chk, chkdetails] = FindCellForExpt(DATA, oldid(1), Clusters,Q,D,I,goodid,state, varargin{:},DATA,'celllist',tC,'nocheck');
                if ~chkdetails.found
                end
                if isfield(chkdetails,'aligndata')
                end
            end
            nbad = 0;
        else
            p = details.nearest(2);
            cl = details.nearest(3);
            if p > 0
                details.xc = details.pxc(p);
                if strcmp(details.fail,'Is')
                    matches(ne).good = -1;
                elseif strcmp(details.fail,'Dr')
                    matches(ne).good = -2;
                elseif strcmp(details.fail,'xc')
                    matches(ne).good = -3;
                else
                    matches(ne).good = 0;
                end
            else
                matches(ne).good = -4;
            end
            nbad = nbad+1;
            if bade == 0
                bade = e;
                Cd{1} = details;
            end
        end
        matches(ne).pos= [e p cl];
        if details.id(2) > 0
            matches(ne).xc = details.xc;
            matches(ne).dropi = D(e,p,cl);
            matches(ne).I = I(e,p,cl);
            matches(ne).xcs = details.probexcs(~isnan(details.probexcs)); %clusters on best probe
            matches(ne).probexc = details.pxc; %max for each probe
        end
        if isfield(details,'aligndata') && ~isempty(details.aligndata)
            matches(ne).aligndata = details.aligndata;
        elseif isfield(chkdetails,'aligndata')
            matches(ne).aligndata = chkdetails.aligndata;
        else
            matches(ne).aligndata = [];
        end
        if length(matches(ne).pos) < 3
            fprintf('short pos\n');
        end
    end
    if nbad > 2 %4 missing expts in a row is time to stop. Can always merge later
        break;
    end
end
goode = id(1);
goodid = tid;
bade = 0;
nbad = 0;
for e = id(1)-1:-1:1;
    ne = ne+1;
    if ismember(e,DATA.autolist.eid)
        chkdetails = [];
        cid = find(C>0);
        tC(cid) = C(cid);
        [x, details] = FindCellForExpt(DATA, e, Clusters,Q,D,I,goodid,state,varargin{:},DATA.ArrayConfig,'celllist',tC);
        if details.newclusters
            Clusters = getappdata(DATA.toplevel,'NewClusters');
            fixdetails.newclusters = fixdetails.newclusters+1;
        end
        if details.fixfit
            fixdetails.fixlist(end+1,:) = [details.id(1:2) details.fixfit 3 now];
            fixdetails.fixfit(details.id(1),details.id(2)) = details.fixfit;            
        end

        oldx = squeeze(CellList(e,:,:));
        [ap,acl] = find(x >0);
        x(abs(oldx)>0) = 0;
        C(e,:,:) = x;
        if details.found
            [p,cl] = find(x >0);
            if isempty(p)
                p = ap(1);
                cl = acl(1);
                matches(ne).good = 2; %duplicate
            else             
            matches(ne).good = 1;
            goode = e;
            goodid = [e p cl cellid];
            end
            nbad = 0;
        else 
            p = details.nearest(2);
            cl = details.nearest(3);
            if p > 0
                details.xc = details.pxc(p);
            else
                details.xc = NaN;
            end
            matches(ne).good =0;
            if bade == 0
                bade = e;
                Cd{2} = details;
            end
            nbad = nbad+1;
        end
        matches(ne).pos= [e p cl];
        if p > 0
            matches(ne).xc = details.xc;
            matches(ne).dropi = D(e,p,cl);
            matches(ne).I = I(e,p,cl);
            matches(ne).xcs = details.probexcs(~isnan(details.probexcs)); %clusters on best probe
            matches(ne).probexc = details.pxc; %max for each probe
        end
        if isfield(details,'aligndata') && ~isempty(details.aligndata)
            matches(ne).aligndata = details.aligndata;
        elseif isfield(chkdetails,'aligndata')
            matches(ne).aligndata = chkdetails.aligndata;
        else
            matches(ne).aligndata = [];
        end
        if length(matches(ne).pos) < 3
            fprintf('short pos\n');
        end
    end
    if nbad > 2 %4 missing expts in a row is time to stop. Can always merge later
        break;
    end
end



function [C, details] = FindCellForExpt(DATA, e, Clusters, Q,D,I, id, state, varargin)

corrth  = [0.98 0.95];
isoth = 2.5;
dropth = 1.5;
cargs = {'noellipse'};
verbose = [1 0];
if isfield(DATA,'ArrayConfig')    
    ArrayConfig = DATA.ArrayConfig;
    cargs = {cargs{:} ArrayConfig};
end
if isfield(DATA.autolist,'corrproberange')
    cargs = {cargs{:} 'proberange' DATA.autolist.corrproberange};
end
if isfield(DATA.autolist,'corrth')
    corrth = DATA.autolist.corrth;
end

checksuspicious = 1;
applyfixes = state.applyfixes;
idetails.newclusters = 0;
idetails.fixfit = 0;
idetails.nearest = [0 0 0];
idetails.id = [0 0 0];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'corrth',6)
        j = j+1;
        if length(varargin{j}) ==1
            corrth(1) = varargin{j};
        else
            corrth = varargin{j};
        end
    elseif strncmpi(varargin{j},'applyfixes',7)
        applyfixes = 1;
    elseif strncmpi(varargin{j},'dropth',6)
        j = j+1;
        dropth = varargin{j};
    elseif strncmpi(varargin{j},'celllist',6)
        j = j+1;
        CellList = varargin{j};
    elseif strncmpi(varargin{j},'isoth',6)
        j = j+1;
        isoth = varargin{j};
    elseif strncmpi(varargin{j},'nocheck',6)
        checksuspicious = 0;
    elseif strncmpi(varargin{j},'silent',6)
        verbose = [0 0];
    elseif strncmpi(varargin{j},'wide',4)
        cargs = {cargs{:} 'delays' 5};
    elseif strncmpi(varargin{j},'verbose',4)
        verbose = [1 1];
    end
    j = j+1;
end

if length(id) > 3
    cellid = sprintf('%d',id(4));
    cellno = id(4);
else
    cellid = '?';
    cellno = 0;
end

C = zeros(size(Q,2),size(Q,3));
    details.found = 0;
    details = CopyFields(details,idetails,'-ifnew');
A = PC.GetClusterInfo(Clusters,id(1),id(2),id(3),'allexpt');
if isempty(A.times) && ~isfield(A,'prpobe')
    fprintf('Cluster E%dP%d.%d empty\n',id);
    return;
end

if e <1 || e > size(Q,1)    
    return;
end
if ~isempty(ArrayConfig)
    probes = array.FindProbe(ArrayConfig,id(2),2);
else
    probes = 1:size(Q,2);
end
    for p = probes(:)' %all probes with 2
        for c = 1:size(Q,3) %all clids
            B = PC.GetClusterInfo(Clusters,e,p,c,'allexpt');
            [xc(c), details] = PC.ShapeCorr(A,B,'drifts',1,cargs{:});
            if isfield(details,'amp')
                amps(c) = details.amp;
            else
                amps(c) = NaN;
            end
            ppx(c) = details.probeshift;
            ppt(c) = details.timeshift;
            alldetails{p,c} = details;
            alldetails{p,c}.xc = xc(c);
        end
        clsxc{p} = xc;
        clsamp{p} = amps;
        [pxc(p), pcl(p)] = max(xc);
        if ppx(pcl(p)) ~= 0
            dp = ppx(pcl(p));
            dt = ppt(pcl(p));
            if pxc(p) > corrth(2);
                fprintf('E%dShape Match for P%d on P%d=%d\n',e,id(2),p,p-dp);
            end
        else
            dp = 0;
            dt = 0;
        end
        timeshift(p) = dt;
        probeshift(p) = dp;
    end
    B.exptno = e;
    details.found = 0;
    [a,b] = max(pxc);
    details.probexcs = clsxc{b};
    details.probeamps = clsamp{b};
    details.pxc = pxc;
    details.ppx = ppx;
    details.isolation = I(e,b,pcl(b));
    details.dropi = D(e,b,pcl(b));
    details.id = [e b pcl(b)];
    details.probeshifts = probeshift; %shift used for matching each probe
    details.aligndata = [];
    details.goodmatch = 0;
    details.xc = a;
    details.newclusters = 0;
    details.probeshift = probeshift(b);
    details.timeshift = timeshift(b);
    
    mx.toplevel = DATA.toplevel;
    mx.cid = [id(3) pcl(b)];
    mx.probeshift = probeshift(b);
    mx.cell = id(4);

%b is probe, pcl(b) is cluster of highest corr match
    pos = [id(1:2) 1; e b 1];
    mx.pos = pos;
    details.nearest = [e b pcl(b)];
    if D(e,b,pcl(b)) > dropth & I(e,b,pcl(b)) > isoth && a > corrth(1)
        details.found = details.found+1;
        if abs(details.probeshift) > 0  % Take care here. Don't switch proibes 
              tpos = pos;
              tpos(2,2) = id(2);%unless sures
              if id(2) == b %spike changed but trigger probe did not
              else
                  Aa = PC.GetClusterInfo(Clusters,tpos,'clst','addfit');
                  [M, x] = clust.MahalMatrix(Aa,'find',id(3));
                  mx = CopyFields(mx,x);
                  mx.probeshift = 0;
                  [details.goodmatch, xx] = CheckMatch(Aa,mx, CellList);
                  if details.goodmatch > 0
                      b = id(2);
                      pcl(b) = mx.bestmatch(1);
                  end
              end
        end
        if sum(details.probexcs > 0.9) > 1 && checksuspicious %two matches
            cid = find(details.probexcs > 0.9);
            amps = abs(details.probeamps-1);
            [~,did] = min(amps(cid));
%checksuspicous == 0 on recursive calls            
            xpos = pos;
            xpos(1,3) = id(3); %will differ if changed probe above
            xpos(2,3) = pcl(b);
            xpos(2,2) = b;
            xpos(1,3) = 1;
            xpos(2,3) = 1;
            Aa = PC.GetClusterInfo(Clusters,xpos,'clst','addfit');
            if state.watchprogress && checksuspicious
                mx.pos = xpos;
                DATA.autolist.CellList = CellList;
                PC.CompareQuickSpks(DATA, xpos,'celllist',CellList);
                Xa = PC.GetClusterInfo(Clusters,xpos,'clst','addfit');
                [M, x] = clust.MahalMatrix(Xa,'plot',[id(3) pcl(b)]);
                drawnow;
            else
                [M, x] = clust.MahalMatrix(Aa);
                mx.pos = xpos;
            end
            mx.M = M;
            mx = CopyFields(mx,x);
            mx.cid =[id(3) pcl(b)];
            [details.goodmatch, xx] = CheckMatch(Aa,mx, CellList);
            details.bettermatch = xx.bettermatch;
            if sum(xx.othermatch > 0) > 2 %two matches
                details.goodmatch = -1; %go to aligh fits below
            end
            if details.goodmatch == -1 || details.bettermatch > 0
                if state.watchprogress
                    clust.MahalMatrix(Aa,'plot',mx.cid);
                end
                if details.bettermatch > 0                    
                    d(1,:) = M(mx.cid(1),mx.cid(2),:);
                    a = find(details.probexcs >0.9);
                    [~,xid] =  sort(details.probexcs(a),'descend');
                    for k = 2:length(xid)
                        d(k,:) = M(mx.cid(1),a(xid(k)),:);
                    end
                    [~,besti] = min(nanmean(d));
                    if diff(mean(d(1:2,:)')) < 0 || sum(diff(d(1:2,:)) <0) > 2 || besti > 1
                        if a(xid(2)) ~= details.bettermatch
                            details.bettermatch = a(xid(2));
                        else %bettermatch is also second highest xcorr
%?? Check that second xc isn't too much lower?                             
                            pcl(b) = details.bettermatch;
                        end
                        pcl(b) = details.bettermatch;
                    else
                        details.bettermatch = -details.bettermatch; %not really true
                    end
                else
                    [~,fix] = clust.AlignAutoFits(Aa,'checkcl',id(3),cargs{:});
                    details.aligndata = fix;
                    details.aligndata.type = 'alignfit';
                    if ismember(fix.usefix,[1 2]) %not ready to use fix 3(ellipses) yet
                        f = fix.newfit(fix.usefix);
                        if f.matchdistance(1) < 4 && applyfixes
                            if fix.usefix == 1
                                cpos = pos(2,1:2);
                            else
                                cpos = pos(1,1:2);
                            end
                            res = PC.CallAllVPcs(DATA.toplevel, cpos(1),cpos(2),'setfit',f.fitnumber);
                            cpfields = clust.autofitfields;
                            Clusters = clust.Copy(Clusters, res.cluster, cpos,cpfields);
                            Clusters{cpos(1)}{cpos(2)}.autofixed = f.fitnumber;
                            setappdata(DATA.toplevel,'NewClusters',Clusters);
                            details.newclusters = 1;
                            details.fixfit = f.fitnumber;
                        end
                    end
                end
            end
            if cid(did) ~= pcl(b) || sum(amps(cid) > 0.2) > 1 && details.newclusters ==0
            Ba = PC.GetClusterInfo(Clusters,[e b 1],'clst','addfit');
            Aa = PC.GetClusterInfo(Clusters,[id(1) id(2) 1],'clst','addfit');
            
            [~,details.aligndata] = clust.AlignAutoFits({Aa Ba}, 'checkcl', id(3),'drifts',1,cargs{:});
            details.aligndata.pos = [id(1:3); e b 1];
            details.aligndata.type = 'alignfit';
            drawnow;
            if details.aligndata.match ~= pcl(b) 
                if probeshift(b) == 0
                    fprintf('May not be best match at E%dP%d,%d\n',e,Aa.probe,Ba.probe);
                    details.aligndata.aligned = 100 * id(3);
                else
                    fprintf('Match Needs Probe Shift. E%dP%d,%d\n',e,Aa.probe,Ba.probe);
                end
            end
            end
        elseif sum(details.pxc > 0.9) > 1 && pxc(id(2)) > 0.9 && b ~= id(2)% match on > 1 probe
            pid = find(details.pxc > 0.9);
            details.altprobe = id(2);
            pos(2,2) = pos(1,2); %check matches on same probe
            Aa = PC.GetClusterInfo(Clusters,pos,'clst','addfit');
            [M, x] = clust.MahalMatrix(Aa);
            mx.M = M;
            mx.pos = pos;
            mx.cid = [id(3) pcl(pos(1,2))]; %best match on same probe
            mx.probeshift = 0;
            [details.goodmatch, xdetails] = CheckMatch(Aa,mx, CellList);
            details.bettermatch = xdetails.bettermatch;
            if sum(xdetails.othermatch > 0) > 2 %two matches
                details.goodmatch = -1; %go to aligh fits below
            end
            if details.goodmatch == -1 || xdetails.bettermatch > 0
                Aa = PC.GetClusterInfo(Clusters,pos,'clst','addfit');
                [M, mdetails] = clust.MahalMatrix(Aa);
            else %match on same probe seems OK
                xm(1,:) = squeeze(M(mx.cid(1),mx.cid(2),:));
                pos(2,2) = b; %best probe for corr
                Aa = PC.GetClusterInfo(Clusters,pos,'clst','addfit');
                [M, xx] = clust.MahalMatrix(Aa);
                xm(2,:) = squeeze(M(mx.cid(1),pcl(b),:));
                xm(isnan(xm)) = 100; %punish NaNs
                dratio = xm(1,:)./xm(2,:);
                dropi(1) = D(pos(2,1),id(2), mx.cid(2)); %best on same probe
                dropi(2) = D(pos(2,1),b, pcl(b));
                isol(1) = I(pos(2,1),id(2), mx.cid(2)); %best on same probe
                isol(2) = I(pos(2,1),b, pcl(b));
                details.aligndata.pos = pos;
                details.aligndata.goodmatch = details.goodmatch;
                details.aligndata.bettermatch = details.bettermatch;
                details.aligndata.mahal = xm;
                details.aligndata.dropi = dropi;
                details.aligndata.isolation = isol;
                details.aligndata.probes = [id(2) b; mx.cid(2) pcl(b)];
                details.aligndata.type = 'mahalmatrix';
                if sum(dratio < 1.5) > 3 && (dropi(1) > 2 || diff(dropi) < 0.5) ... %seems good enough
                        && (isol(1) > isol(2) * 0.8 || isol(1) > 4)
                    ob = b;
                    details.altprobe = ob;
                    b = id(2);
                    pcl(b) = mx.cid(2);
                    details.xc = pxc(b);
                    fprintf('Sticking with %d: (%.3f:I%.2f,D%.2f) rather than %d (%.3f:I%.2f,D%.2f)\n',...
                    b,pxc(b),isol(1),dropi(1),ob,pxc(ob),isol(2),dropi(2));
                    
                end
                details.aligndata.probe = b;
            end
            if probeshift(b) ~= 0 %also seems that spike has really moved
            end
        end
        C(b,pcl(b)) = 1;
%Set these again in case b or pcl(b) was changed        
        details.isolation = I(e,b,pcl(b));
        details.dropi = D(e,b,pcl(b));
        details.id = [e b pcl(b)];
        details.probeshift = probeshift(b);
        details.timeshift = timeshift(b);
        if verbose(1)
            fprintf('Expt%dP%d.%d  matches cell%s on E%dP%d.%d. Xcorr %.4f\n ',e,b,pcl(b),cellid,id(1),id(2),id(3),a);
        end
    else
        if verbose (2)
            fprintf('Expt%d no match for cell%s on E%dP%d.%d. Closest%.3f: P%dcl%d: ',e,cellid,id(1),id(2),id(3),a,b,pcl(b));
        if D(e,b,pcl(b)) <= dropth
            fprintf(' dropi %.2f\n',D(e,b,pcl(b)));
            details.fail = 'Dr';
        elseif I(e,b,pcl(b)) <= isoth
            fprintf(' isolation %.2f\n',I(e,b,pcl(b)));
            details.fail = 'Is';
        elseif a <= corrth(1)
            details.fail = 'xc';
            fprintf(' ShapeCorr %.4f\n',a);
        else
            fprintf('??\n');
            details.fail = '??';
        end
        end
    end
if checksuspicious && details.found %only check if think its a match
    clear B;
    if details.found && details.id(2) ~= id(2) %change of probe
%Check Match in new expt(new probe) with orignal probe new expt
%high efficacy means reall is the same
        B = PC.GetClusterInfo(Clusters,[details.id; e id(2) pcl(id(2))]);
        if details.pxc(id(2)) > corrth(2)
            eff = CalcEfficacy(B{1},B{2});
            eff = min(eff(5:6));
            if details.pxc(id(2)) > corrth(1) || details.pxc(id(2))> pxc(b).*0.98
                fprintf('E%d->%d Cell Changing Probe(%d) despite good corr on P%d of %.3f but eff%.3f\n',id(1),e,details.id(2),id(2),details.pxc(id(2)),eff);
                if state.watchprogress
                    xpos = [id(1:3); e id(2) pcl(id(2))];
                    PC.CompareQuickSpks(DATA, xpos,'celllist',CellList);
                end
                if pxc(id(2))./pxc(b) > 0.98 %similar scores also - stick with original
                    fprintf('Staying with Original Probe\n');
                    b = id(2);
                    details.isolation = I(e,b,pcl(b));
                    details.dropi = D(e,b,pcl(b));
                    details.id = [e b pcl(b)];
                    details.probeshift = probeshift(b);
                    details.timeshift = timeshift(b);
                    details.xc = pxc(b);
                end
            else
                fprintf('E%d->%d Cell Changing Probe despite moderate corr on P%d of %.3f but eff%.3f\n',id(1),e,id(2),details.pxc(id(2)),eff);
            end
        end
    else
        B{1} = PC.GetClusterInfo(Clusters,[details.id]);
    end
    fprintf('Checking match other way: E%dP%d.%d -> E%dP%d.%d\n',details.id(1:3),id(1:3));
%always check that teh match works both ways    
    [~, chkdetails] = FindCellForExpt(DATA, id(1), Clusters, Q,D,I, [details.id id(4)], state, varargin{:},'nocheck','silent','celllist',CellList,'ischeck');
    if sum(chkdetails.id(1:3)) == 0 
        fprintf('E%DP%D.%d produces bad Check',details.id(1:3));
    elseif sum(chkdetails.id(1:3) == id(1:3)) < 3
        Ca = PC.GetClusterInfo(Clusters,chkdetails.id);
        [eff, ex] = CalcEfficacy(A,Ca);
        %if new cell matches two cells in previous expt, those two may be partial duplicates
        %if efficacy between them is strong but asymmetrical, then this cell is not
        %a match for both. The correct match should pass the Expt1->Expt2->Expt1 loop
        %that this failed.
        fprintf('%s:',datestr(now));
        if isempty(Ca.times)
            fprintf('E%DP%D.%d no longer defined',chkdetails.id(1:3));
        elseif min(eff(5:6)) < 0.2 & max(eff(5:6)) > 0.5
            details.match = 0;
            details.fail = 'Dsplit';
            fprintf('!!!!E%dP%d.%d(%d)->E%d.P%d.%d(%d)->E%dP%d.%d(%d) Eff%.2f and%.2f!!!!\n',id,A.ncut,details.id,B{1}.ncut,chkdetails.id,Ca.ncut,eff(5:6));
        else
            fprintf('Duplicate match OK: E%dP%d.%d(%d)->E%d.P%d.%d(%d)->E%dP%d.%d(%d) Eff%.2f and%.2f\n',id,A.ncut,details.id,B{1}.ncut,chkdetails.id,Ca.ncut,eff(5:6));
        end
    end
    if isempty(details.aligndata) && isfield(chkdetails,'aligndata') 
        details.aligndata = chkdetails.aligndata;
    end
end
details = CopyFields(details,idetails,'-ifnew');
    
    
function [goodmatch, details] = CheckMatch(C,x, CellList,varargin)
    
 goodmatch = -1;
 details.bettermatch = 0;
 xmatch = x.M(x.cid(1),x.cid(2),:);
 details = CopyFields(details,x,{'pos' 'probeshift'});
 
 cls = unique(C{2}.clst);
 details.othermatch = [];
 if ~cellstrcmp('one',varargin)
     tmpx = x;
     for j = 2:length(cls)
         tmpx.cid(2) = cls(j)-1;
         details.othermatch(j-1) = CheckMatch(C,tmpx,CellList,'one');
     end
 end
 
 if x.probeshift ~= 0
     pos = x.pos;
%See if matches on unshifted probes are duplicate cells
     pos(2,2) = pos(1,2);
     Aa = PC.GetClusterInfo(x.toplevel,pos);
     [M, mdetails] = clust.MahalMatrix(Aa,'shift');
     minM = min(M,[],3);
     [d(1),match(1)] = min(minM(x.cid(1),:));
     
     pos = x.pos;
     pos(1,2) = pos(2,2);
     Aa = PC.GetClusterInfo(x.toplevel,pos);
     [M, mdetails] = clust.MahalMatrix(Aa);
     minM = min(M,[],3);
     [d(2),match(2)] = min(minM(:,x.cid(2)));
     if abs(CellList(x.pos(1,1),x.pos(2,2),match(2))) == x.cell
         goodmatch = 2;
     elseif abs(CellList(x.pos(2,1),x.pos(1,2),match(1))) == x.cell
         goodmatch = 1;
     end
     if goodmatch > 0
         fitstr{1} = sprintf('Good Match - probe drift');
     end
 elseif xmatch(2) < 4.5 && xmatch(3) < 4.5
     goodmatch = 3;
     fitstr{1} = sprintf('Good Match %s',sprintf(' %.2f',xmatch));
 elseif xmatch(1) < 4.5 && xmatch(4) < 4.5
     fitstr{1} = sprintf('Good Match %s',sprintf(' %.2f',xmatch));
     goodmatch = 4;
 elseif sum(xmatch <4.5) > 1
     goodmatch = 5;
 end
 minM = squeeze(min(x.M,[],3));
 maxM = squeeze(max(x.M,[],3));
 [a, ai] = min(minM(x.cid(1),:));
 [b, bi] = min(maxM(x.cid(1),:));
 if ai == bi && ai == x.cid(2) && goodmatch == 0%this is best match
     goodmatch = 6;
 end
 details.bestmatches = [ai bi];
 
 if ai ~= x.cid(2) && bi~= x.cid(2)
     if a < b
         details.bettermatch = ai;
     else
         details.bettermatch = bi;
     end
 end
 
 
    function C = FindMuCells(DATA, CellList, Clusters, Q, D,I, state, varargin)
        cells = unique(CellList(CellList>0));
        muCellList = zeros(size(CellList));
        for j = 1:length(cells)
            cid = find(CellList == cells(j));
            [e,p,cl] = ind2sub(size(CellList),cid);
            gaps = setdiff(e(1):e(end),e);
            for k = 1:length(gaps)
                a = find(e < gaps(k));
                a = a(end);
                tid = [e(a) p(a) cl(a) cells(j)];
                [x, details] = FindCellForExpt(DATA,gaps(k), Clusters,Q,D,I,tid,state,'wide','corrth',[0.95 0.9],'dropth',1,'isoth',2,DATA,'celllist',CellList);
                xcell = CellList(details.id(1),details.id(2),details.id(3));
                if details.found && xcell == 0
                    muCellList(details.id(1),details.id(2),details.id(3)) = cells(j);
                else
                    [mx,mp] = max(details.pxc);
                    fprintf('E%dP%dC%d Cell %d best Corr %.2f with E%dP%dC%d cell %d',...
                        e(a),p(a),cl(a),cells(j),mx,details.id,xcell);
                    if mx > 0.95
                        muCellList(details.id(1),details.id(2),details.id(3)) = cells(j);
                        if details.dropi > 1.2 && details.isolation > 2.1
                            fprintf('Added to MUlist\n');
                        else
                            fprintf('Dropi %.2f, Isolation %.2f\n',details.dropi,details.isolation);
                        end
                    else                        
                        PC.CompareSpikeShape(Clusters, [e(a) p(a) cl(a); details.id],'tag','CompareShape');
                        fprintf('NonMatch\n');
                    end
                end
            end
        end
        C = muCellList;
        
    function [C, details] = CheckDuplicates(DATA, CellList, D, I, Clusters, varargin)
        
        details = [];
        cellgroups = {};
        G = getappdata(DATA.toplevel,'AllDuplicates');
        C = CellList;
        ncells = max(C(:));
        ng = 0;
        for j = 1:length(G)
            if isfield(G{j},'nduplicates') && G{j}.nduplicates > 0
                for k = 1:length(G{j}.groups)
                    if length(G{j}.groups{k}) > 1
                        ng = ng+1;
                        groups{ng} = G{j}.groups{k};
                        probes{ng} = floor(G{j}.probes(G{j}.groups{k}));
                        dprobes{ng} = G{j}.probes(G{j}.groups{k});
                        clid{ng} = round(rem(G{j}.probes(G{j}.groups{k}),1).*10);
                        gid(ng) = j;
                        kid(ng) = k;
                    end
                end
            end
        end
        for j = 1:ng
            expts = ones(size(probes{j})) .* gid(j);
            id = sub2ind(size(CellList),expts, probes{j}, clid{j});
            cells{j} = CellList(id);
            cid = find(cells{j} > 0);
            c = cells{j}(cid);
            xp = cat(1,G{gid(j)}.xcorrs.probe);
            if length(c) > 1
                a = find(sum(ismember(xp,dprobes{j}(cid)),2)==2 & diff(xp,[],2) ~= 0);
                xeff = cat(1,G{gid(j)}.xcorrs(a).efficacy);
                if max(min(xeff(:,5:6))) < 0.2
                    fprintf('E%d still has Partial duplicate cells%',gid(j));
                else
                    fprintf('E%d still has duplicate cells%',gid(j));
                end
                    
                for k = 1:length(c)
                    fprintf(' %d (%.1f)',c(k),dprobes{j}(cid(k)));
                end
                fprintf('\n');
                for k = 1:length(a)
                    x = G{gid(j)}.xcorrs(a(k));
                    fprintf('%.1f <-> %.1f %.2f,%.2f\n',x.probe(1),x.probe(2),x.efficacy(5:6));
                end
            end
        end
        
function [CellList, duplicates] = OldFindDuplicates(CellList, I, D)
    duplicates = [];
    for e = 1:size(CellList,1);
        C = squeeze(CellList(e,:,:));
        [p,c] = find(C>0);
        E = [];
        for j = 1:length(p)
            A = PC.GetClusterInfo(Clusters,[e p(j) c(j)]);
            for k = 1:j-1
                B = PC.GetClusterInfo(Clusters,[e p(k) c(k)]);
                eff = CalcEfficacy(A,B);
                E(j,k) = eff(1);
                Eb(j,k) = eff(2);
            end
        end
        if length(p) > 1
            [a,b] = find(E > 0.2);
            for j = 1:length(a)
                cella =  [CellList(e,p(a(j)),c(a(j))) I(e,p(a(j)),c(a(j))) D(e,p(a(j)),c(a(j)))];
                cellb =  [CellList(e,p(b(j)),c(b(j))) I(e,p(b(j)),c(b(j))) D(e,p(b(j)),c(b(j)))];
                eff = [E(a(j),b(j)) Eb(a(j),b(j))];
                if cella(2) > cellb(2) && cella(3) > 2
                    fprintf('E%dCell%d (%.2f,%.2f) is copy of cell%d (%.2f,%.2f) eff %.3f,%.3f\n',e,cella,cellb,eff);
                    CellList(e,p(b(j)),c(b(j))) = -CellList(e,p(a(j)),c(a(j)));
                    duplicates(end+1,:) = [e cella(1) cellb(1)];
                elseif cellb(3) > 2 && cellb(2) > 3
                    fprintf('E%dCell%d (%.2f,%.2f) is copy of cell%d (%.2f,%.2f) eff %.3f,%.3f\n',e,cellb,cella,eff);
                    CellList(e,p(a(j)),c(a(j))) = -CellList(e,p(b(j)),c(b(j)));
                    duplicates(end+1,:) = [e cellb(1) cella(1)];
                else
                    fprintf('Cell%d ,%d duplicates\n',cellb,cella);
                end
            end
        end
    end
    
    
    
    
function [hcount, dup, xData] = CountDuplicates(CellList, G, Clusters, varargin)    
    hcount = [];
    dup = [];
    interactive = 0;
strargs = {};    
checkoverlap = 0;
xData = [];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'checkoverlap',6)
        checkoverlap = 1;
        interactive = 1;
    elseif ischar(varargin{j})
        strargs = {strargs{:} varargin{j}}; 
    end
    j = j+1;
end
        ncells = max(CellList(:));
        ng = 0;
        for j = 1:length(G)
            if isfield(G{j},'nduplicates') && G{j}.nduplicates > 0
                for k = 1:length(G{j}.groups)
                    if length(G{j}.groups{k}) > 1
                        ng = ng+1;
                        groups{ng} = G{j}.groups{k};
                        cellgroups{ng} = G{j}.cellid(G{j}.groups{k});
%use existing correlatinos, but allow for changes in celllist numbers                        
                        probes{ng} = floor(G{j}.probes(G{j}.groups{k}));
                        clid{ng} = round(rem(G{j}.probes(G{j}.groups{k}),1).*10);
                        gid(ng) = j;
                        expts = ones(size(probes{ng})) .* j;
                        id = sub2ind(size(CellList),expts, probes{ng}, clid{ng});
                        c = CellList(id);
                        c = c(c>0);
                        cellgroups{ng} = c;
                    end
                end
            end
        end
        hcount = zeros(ncells);
        pcount = hcount;
%cellgroups are list of cells that duplicated each other

        for j = 1:length(cellgroups)
            e = G{gid(j)}.expt;
            e = gid(j);
            for k = 1:length(cellgroups{j});
                for t = 1:k-1
                    a = cellgroups{j}(k);
                    b = cellgroups{j}(t);
                    if b == a %self duplicate. Error, but should not count as dup
                        fprintf('!!!!E%d Cell%d is Self Duplicate\n',e,a);
                    else
                        isdup = PC.IsSameCell(PC.FindXCorr(G{gid(j)}.xcorrs,groups{j}([k t])));
                        if isdup ==2 && sum(strcmp('partial',strargs));
                            dup(e,a) = b;
                            dup(e,b) = a;
                        else
                            dup(e,a) = b;
                            dup(e,b) = a;
                        end
                        exptgroup(e) = gid(j);
                        if isdup == 1
                            if b > a
                                hcount(b,a) = hcount(b,a)+1;
                            else
                                hcount(a,b) = hcount(a,b)+1;
                            end
                        else
                            if b > a
                                pcount(b,a) = pcount(b,a)+1;
                            else
                                pcount(a,b) = pcount(a,b)+1;
                            end
                        end
                    end
                end
            end
        end
%For cells that are clear duplicates in some places, count partial duplicates too        
if sum(strcmp('partial',strargs))
    hcount = pcount;
else
    id = find(hcount >1);
    hcount(id) = hcount(id)+pcount(id);
end
%for cells that only show partial duplication, surely one should go. ? deal
%with these cases first, or last? 
%?mark for resolution by user? 
   
    if checkoverlap
       [a,b] = find(hcount > 0);
       allxc = [];
       alleff = [];
       allxid = [];
       allcid = [];
       npair = 0;
       np = 0;
       for j = 1:length(a)
           if a(j) == b(j)
           else
               [ea, pa, cla] = celllist.find(CellList, a(j));
               [eb, pb, clb] = celllist.find(CellList, b(j));
               [aid, ia, ib] = intersect(ea,eb); %expts with both
               id = find(dup(:,a(j)) > 0 & dup(:,b(j)) >0); %expts with dup
               nondupid = setdiff(aid,id);
               dupid = id;
               if isempty(id)
                   fprintf('!!!%d->%d(%d) empty list!!!!\n');
               end
               isdup = [];
               for k = 1:length(aid)
                   e = aid(k);
                   ca = pa(ia(k)) + cla(ia(k))/10;
                   cb = pb(ib(k)) + clb(ib(k))/10;
                   p = cat(1,G{e}.xcorrs.probe);                   
                   xid = find((p(:,1) == cb & p(:,2) == ca) | (p(:,1) == ca & p(:,2) == cb));
                   [isdup(k) score(k,:)]= PC.IsSameCell(G{e}.xcorrs(xid));
               end
               if ~isempty(nondupid) && sum(isdup ==1) == 0
%all teh duplicates are partil.  Certainly don't want automatic process to treat as duplicates                   
%?look at the low scores?
                  fprintf('Cell%d partial duplicate of %d in %d expts, but no dup in %d others\n',...
                        a(j), b(j), length(id),length(nondupid));
                  hcount(a(j),b(j)) = 0;
               elseif ~isempty(nondupid) 
                   np = np+1;
                    fprintf('Cell%d duplicates %d in %d expts, but NOT in %d others\n',...
                        a(j), b(j), length(id),length(nondupid));
                    if interactive
                        allid = cat(1,id, nondupid);
                        for k = 1:length(allid);
                            npair = npair+1;
                            e = allid(k);
                            ca = find(G{e}.cellid == a(j));
                            cb = find(G{e}.cellid == b(j));
                            cella = G{e}.cells(ca);
                            cellb = G{e}.cells(cb);
                            if abs(cella.p - cellb.p) > 15
                                fprintf('Very Distant Pair\n');
                            end
                            allcid(npair,:) = [cella.eid cella.p cella.cl cellb.eid cellb.p cellb.cl];
                            C = PC.GetClusterInfo(Clusters,G{e}.cells([ca cb]));
                            p = cat(1,G{e}.xcorrs.p);
                            xid = find((p(:,1) == cb & p(:,2) == ca) | (p(:,1) == ca & p(:,2) == cb));
                            alleff(npair) = max(G{e}.xcorrs(xid).efficacy(5:6)); 
                            allxc(npair) = PC.ShapeCorr(C{1},C{2});
                            isdup = ismember(allid(k),nondupid);
                            pairlist(npair,:) = [np isdup a(j) b(j)];
                        end
                    end
                    id = find(pairlist(:,1) == np); %current set
                    xData.range(np,:) = minmax(alleff(id));
                    xData.clusterid = allcid;
                    xData.pairlist = pairlist;
                    xData.alleff = alleff;
                    xData.allxc = allxc;
                    if max(alleff(id)) < 0.5
                        dup(allid,a(j));
                    elseif min(alleff(id)) > 0.5
                        dup(allid,a(j));
                    end
               end
           end
       end
       if interactive && ~isempty(allxc);
           F = GetFigure('Duplicates');
           hold off;
           plot(allxc,alleff,'o','buttondownfcn',{@HitEffScatter});
           xlabel('Shape Corr');
           ylabel('efficacy');
           title('Cell pairs that are duplicate in some expts');
           xData.xy = cat(1,allxc, alleff);
           xData.PCtoplevel = 1; %kludge for now
           set(F,'UserData', xData);
       end
   end
   
function AddEffMenu(F,X)

    it = findobj(allchild(F),'flat','tag','DuplicateMenu');
    if isempty(it)
        it = uimenu(F,'Label','Expts','tag','DuplicateMenu');
    else
        delete(allchild(it));
    end
    for j = 1:size(X,1);
        uimenu(it,'label',sprintf('E%d',X(j,1)),'callback',{@SetEffExpt X(j,1)});
    end
    

function SetEffExpt(src, event,eid)
X = GetDataFromFig(src);
pid= getappdata(gcf,'currentpairs');
id = find(X.clusterid(pid,1) == eid);
e = pid(id);
PlotCellComparison(X, e);
    
        
function HitEffScatter(src, event)

X = GetDataFromFig(src);
[exy, id] = FindNearestPoint(src);
idlist = get(src,'Userdata');
if length(idlist) >= id
    id = idlist(id);
end
pid = find(X.pairlist(:,1) == X.pairlist(id,1)); %other expts with this pair
qid = find(X.pairlist(pid,2) ==1); %nondup
qid = pid(qid);
did = find(X.pairlist(pid,2) ==0); %dup
did = pid(did);
hold off;
aid = setdiff(1:size(X.pairlist,1),pid);
h = plot(X.xy(1,aid,1),X.xy(2,aid),'o','buttondownfcn',{@HitEffScatter});       
set(h,'UserData',aid);
hold on;
h = plot(X.xy(1,qid,1),X.xy(2,qid),'o','buttondownfcn',{@HitEffScatter});       
set(h,'UserData',qid);
h = plot(X.xy(1,did,1),X.xy(2,did),'s','buttondownfcn',{@HitEffScatter});       
set(h,'UserData',did);

G = getappdata(X.PCtoplevel,'AllDuplicates');
[~,sid] = sort(X.clusterid(pid,1));
pid = pid(sid);

for j = 1:length(pid)
    if pid(j) == id %selected point
        fprintf('*');
    end
    PrintDup(X.clusterid(pid(j),:),G);
end
setappdata(gcf,'currentpairs',pid);
AddEffMenu(gcf, X.clusterid(pid,:));
PlotCellComparison(X, id);


function PrintDup(clusterid,G)
    
    e = clusterid(1);
    p = clusterid([2 5]);
    cl = clusterid([3 6]);
    a = find([G{e}.cells.p] == p(1) & [G{e}.cells.cl] == cl(1));
    cella = G{e}.cellid(a);
    b = find([G{e}.cells.p] == p(2) & [G{e}.cells.cl] == cl(2));
    cellb = G{e}.cellid(b);
    plist = cat(1,G{e}.xcorrs.p);
    c = find(plist(:,1) == max([a b]) & plist(:,2) == min([a b]));
    if isempty(c)
        fprintf('E%d Cells %d<->%d\n',e,cella,cellb);
    else
        xc = G{e}.xcorrs(c);
        [isdup,score] = PC.IsSameCell(xc);
        fprintf('E%d Cells %d(%d)<->%d(%d) DUP%d Eff%s\n',e,cella,xc.nspk(1),cellb,xc.nspk(2),isdup,sprintf('%.2f ',xc.efficacy));
    end

function PlotCellComparison(X,id)
    Clusters = getappdata(X.PCtoplevel,'AutoClusters');
    C = PC.GetClusterInfo(Clusters,[X.clusterid(id,1:3); X.clusterid(id,4:6)]);
    exid = X.clusterid(id,[1 4]);
    [xc, details] = clust.xcorr(C);
    [eff, effd] = clust.Efficacy(C);
    GetFigure('xcorr');
    plot(details.xpts, xc);
    title(sprintf('Eff%.2f,%.2f:%.2f,%.2f',details.efficacy,eff(1:2)));
    GetFigure('CompareShape');
    [sxc, details] = PC.ShapeCorr(C, 'delays', 5); %allow change in time, but not probe
    fprintf('Delay %.0f\n',details.timeshift);
    sxc = PC.CompareSpikeShape(C,'delay',[0 details.timeshift]);
    DATA = get(X.PCtoplevel,'UserData');
    DATA = PC.SetValue(DATA,'Expt',exid(1));
    e = clust.Efficacy(C);
    isdup = PC.IsSameCell(e);
    DATA.currentpoint(2) = C{1}.probe(1);
    SetData(DATA);
    
    PC.SetFigure(DATA,DATA.tag.spikes);
    h = PC.QuickSpikes(DATA, [exid(1) C{1}.probe(1)],'autolist');
    title(sprintf('Cell%d E%dP%d.%d %d spikes',X.pairlist(id,3),exid(1),C{1}.probe(1),C{1}.cluster,C{1}.ncut));
    
    PC.SetFigure(DATA,['Compare' DATA.tag.spikes]);
    h = PC.QuickSpikes(DATA, [exid(2) C{2}.probe(1)],'autolist');
    title(sprintf('Cell%d E%dP%d.%d %d spikes',X.pairlist(id,4),exid(2),C{2}.probe(1),C{2}.cluster,C{2}.ncut));
    G = getappdata(X.PCtoplevel,'AllDuplicates');
%    PrintDup(X.clusterid(id,:),G);
%sxc = PC.CompareSpikeShape(C);
        
   
    function [C, details] = FixDuplicateSquares(DATA, CellList, Clusters,Fx)
        
        details.dup = [];
        expts = DATA.autolist.eid;
        for e= expts(:)'
            eC = squeeze(CellList(e,:,:));
            [p,cl] = find(eC ==0);
            [cp, ccl] = find(eC>0);
            for j = 1:length(p)
                for k = 1:length(cp)
                    if Fx.D(e,p(j),cl(j)) > 2.5
                        C =  PC.GetClusterInfo(Clusters,[e p(j) cl(j); e cp(k) ccl(k)]);
                        x = CalcEfficacy(C{1},C{2});
                        if max(x(1:2)) > 0.2
                            cell = CellList(e,cp(k),ccl(k));
                            CellList(e,p(j),cl(j)) = -cell;
                            fprintf('E%dP%d.%d is duplicate of cell %d\n',e,p(j),cl(j),cell);
                            details.dup(end+1,:) = [e p(j) cl(j) cell];
                        end
                    end
                end
            end
          
        end
C = CellList;

        
function [C, details] = CheckOrphanMerge(DATA, C, Clusters, goodcell, badcell)
%Check if goodcell and badcell look like the same based on shape;
%applied after fixing duplicates, so goodcell and badcell should be on
%mutually exclusive experiements
    
    checkpoints = [];
    id = find(C == badcell);
    [e,p,cl] = ind2sub(size(C),id);
    
    gid = find(C == goodcell);
    [ge,gp,gcl] = ind2sub(size(C),gid);
    
    [ebid, epid, ecid] = find(C == -badcell);
    oid = find(ismember(epid,e));
    epid = epid(oid);
    ebid = ebid(oid);
    ecid = ecid(oid);
    %                    oid = find(ismember(ebid,e)); %expts where there are duplicates of badcell
    
    ediff = bsxfun(@minus,ge,e'); %size ge x e
    dcrit = min(abs(ediff(:)))+1;
    dcrit = max([dcrit 3]);
    %ga is a list of all expts containing cell g, where the nearest cell b is < dcrit expts away
    % will get mulitple entries is cell be is 1 expt away and two expts away,
    % but that is OK
    %
    [ga,ba]  = find(abs(ediff) < dcrit & abs(ediff) > 0);
    details.ga = ge(unique(ga));
    ga = unique(ga);
    [oa,ob] = find(ediff ==0); %cells on same expt -
    xcs = [];
    for j = 1:length(ga)
        A = PC.GetClusterInfo(Clusters,ge(ga(j)),gp(ga(j)),gcl(ga(j)),'allexpt');
        B = PC.GetClusterInfo(Clusters,e(ba(j)),p(ba(j)),cl(ba(j)),'allexpt');
        [xcs(j), xcdetails{1}] = PC.ShapeCorr(A,B,'delays',5,'proberange',5,'drifts',0);
    end
    if isempty(xcs) %No Nearby Expts
        if ~isempty(id)
            fprintf('Cell %d Still %d Expts, not close to Cell %d\n',badcell, length(id),goodcell);
        end
    elseif max(xcs) > 0.95 && min(xcs) >= 0.9 %all expts in ga are good
        oid = sub2ind(size(C),e(ba),p(ba),cl(ba)); %adjacent expts with match
        oid = unique(oid);
        fprintf('Setting Cell %d to %d in %d non-overlapping expts\n',badcell,goodcell,length(oid));
        dup(e(ba),badcell) = 0;
        C(oid) = goodcell;
        C(id(ob)) = -goodcell;
        oid = sum(C(:) == badcell);
        if oid
            fprintf('Cell%d still present in %d other expts',badcell,oid);
        end
        %                       for j = 1:length(oid)
        %                          C(ebid(oid(j)),epid(oid(j)),ecid(oid(j))) = -goodcell;
        %                      end
    elseif max(xcs) > 0.95
        oid = sub2ind(size(C),e(ba),p(ba),cl(ba)); %adjacent expts with match
        checkpoints = [checkpoints oid(xcs < 0.9)'];
        dup(e(ba),badcell) = 0;
        C(oid) = goodcell;
        fprintf('Lowest Corr(%.2f), but Setting Cell %d  matches %d in %d/%d non-overlapping expts\n',min(xcs),badcell,goodcell,sum(xcs>0.95),length(oid));
        C(id(ob)) = -goodcell;
        %                        for j = 1:length(oid)
        %                           C(ebid(oid(j)),epid(oid(j)),ecid(oid(j))) = -goodcell;
        %                        end
        [e,p,cl] = celllist.find(C,badcell);
        ne = length(e);
        if ne
            fprintf('Cell%d still present in %d other expts E%d-%d,P%.1f\n',badcell,ne,min(e),max(e),mean(p));
        end
    else
            fprintf('Cell%d in Expts %s Not a duplicate of %d (best xc %.2f)\n',badcell,sprintf(' %d',e),goodcell,max(xcs));        
    end
    details.checkpoints = checkpoints;
        
function CheckDupChange(xdata, DATA, Clusters, G)
%for pairs where they are duplicate in some expts and not others
%look closely at what happens at the Boudary....

    a = unique(xdata.pairlist(:,3));
    for j = 1:length(a)
        aid = find(xdata.pairlist(:,3) == a(j));
        bcells = unique(xdata.pairlist(aid,4));
        for k = 1:length(bcells)
            dup = [];
            id =find(xdata.pairlist(:,3) == a(j) & xdata.pairlist(:,4) == bcells(k));
            did = find(abs(diff(xdata.pairlist(id,2)))>0); %changes in dup status
            for e = 1:length(did)
                ex(1) = xdata.clusterid(id(did(e)),1);
                ex(2) = xdata.clusterid(id(did(e))+1,1);
                A = PC.GetClusterInfo(Clusters,xdata.clusterid(id(did(e))+[0:1],1:3));
                xc(1) = PC.ShapeCorr(A);
                B = PC.GetClusterInfo(Clusters,xdata.clusterid(id(did(e))+[0:1],4:6));
                xc(2) = PC.ShapeCorr(B);
                if xc(1) < 0.95
                    PC.CompareQuickSpks(DATA,A);
                end
                if xc(2) < 0.95
                    PC.CompareQuickSpks(DATA,B);
                end
                    %                xc = PC.FindXCorr(G{e},
%                dup(e) = PC.IsSameCell(xdata.alleff(id(e),:));
            end
        end
    end
    
    function [C, details] = FixDuplicates(DATA, CellList, D, I, Clusters, varargin)
        
        checkpoints = [];
        details.replaced = [];
        cellgroups = {};
        G = getappdata(DATA.toplevel,'AllDuplicates');
        C = CellList;

        details.initialpartiallist = FindPartialDuplicates(C, G, Clusters);

        [hcount, dup, xData] = CountDuplicates(C ,G,Clusters,'checkoverlap');
        replaced = [];
        replacelist = [];
        [a,b] = find(hcount > 0);
%Hard to keep track of cells as they get replaced because
%original duplicate listting (G) is based on groupings within
%each expt, not just pairised corrs.
%Could go through G each time a cell is replaced and find all
%matches.  But thats a pain.
       details.duplist = hcount;
       lastdup = [];
       while ~isempty(a)
           replaced = [];
           nc = 0;
           for j = 1:length(a)
%if reached a cell htat has already been touched, rebuild the list first
            errs = 0; 
            if sum(ismember([a(j) b(j)],replaced(:)))
                   break;
               end
            if b(j) == a(j)
            else
                nc = nc+1;
                acell = a(j);
                bcell = b(j);
%Get list of expts where these are duplicates of each other.
%N.B. There may be expts where both appear but are not duplicates.
%These will cause troube if a Merge is attempted.
                ea = celllist.find(C, acell);
                eb = celllist.find(C, bcell);
                aid = intersect(ea,eb); %expts with both
                id = find(dup(:,acell) > 0 & dup(:,bcell) >0); %expts with dup
                nondupid = setdiff(aid,id);
                if isempty(id)
                        fprintf('!!!%d->%d(%d) empty list!!!!\n');
                        errs(end+1) = 1;
                end
                nondupid = setdiff(aid,id);
                if ~isempty(nondupid)
                    fprintf('Cell%d duplicates %d in %d expts, but NOT in %d others\n',...
                        acell, bcell, length(id),length(nondupid));
                        errs(end+1) = 2;
                end
                na = sum(C(:) == acell);
                nb = sum(C(:) == bcell);
                clear good;
                clear bad;
                if na < nb
                    good(1:length(id)) = bcell;
                    bad(1:length(id)) = acell;
                    goodcell = bcell;
                    badcell = acell;
                else
                    bad(1:length(id)) = bcell;
                    good(1:length(id)) = acell;
                    goodcell = acell;
                    badcell = bcell;
                end
                    
                fprintf('Cell%d in %dex, Cell%d %dex, %d/%d duplicate. Using cell %d Ex%s\n',...
                a(j),na,b(j),nb,length(id),hcount(a(j),b(j)),goodcell,sprintf(' %d',id));  
                for k = 1:length(id)
                    e = id(k);
                    eC = squeeze(C(id(k),:,:));
                    gid = find(eC == good(k));
                    bid = find(eC == badcell);
                    bid = setdiff(bid,gid); %don't self replace
%A cell might be missing if it was a duplicate with someone
%else and so was reset in a previous loop
                    if ~isempty(bid) && ~isempty(gid)
                        if length(gid) > 1
                            fprintf('Duplicate Good cell %d in Expt %d\n',good(k),e)
                            ind2sub(size(eC),gid);
                            errs(end+1) = 3;
                        end
                        if length(bid) > 1
                            fprintf('Duplicate Bad cell %d in Expt %d\n',bad(k),e)
                            ind2sub(size(eC),gid);
                            errs(end+1) = 4;
                        end
                        [p, cl] = ind2sub(size(eC),gid(1));
                        [p(2:length(bid)+1), cl(2:length(bid)+1)] = ind2sub(size(eC),bid);
                        best(k) = good(k);
                        if D(e,p(2),cl(2)) > D(e,p(1),cl(1)) && I(e,p(2),cl(2)) > I(e,p(1),cl(1))
                            if I(e,p(1),cl(1)) < 4 || D(e,p(1),cl(1)) < 2
                                best(k) = bad(k);
                                bad(k) = good(k);
                                good(k) = best(k);
                                p = fliplr(p);
                                cl = fliplr(cl);
                            end
                        end
                        
%p(2) cl(2) now slated for replacement.  Check that it was not considered the best cell
%in GroupDuplicates                        
                        probe = p + cl./10;
                        ca = find([G{e}.cells.p] == p(1) & [G{e}.cells.cl] == cl(1));
                        cb = find([G{e}.cells.p] == p(2) & [G{e}.cells.cl] == cl(2));
                        cid = find(sum(cellmember([ca cb],G{e}.groups)) > 1); %should be 2
                        if ~isempty(cid) && cb == G{e}.bestcell(cid);
                            fprintf('FixDup:  Bad cell %d was best in Group Expt %d\n',bad(k),e)                            
                            errs(end+1) = 5;
                        end
%anthing currrently recorded as a copy of badcell should be set to a copy of goodcell                        
                        [bpid, bcid] = find(squeeze(C(id(k),:,:)) == -badcell);
                        for c = 1:length(bpid)
                            C(id(k),bpid(c),bcid(c)) = -abs(goodcell);
                        end
                        replaced(end+1,:) = [badcell goodcell];
%replacelist records original cell# in locattions that get marked as duplicates                        
                        rid = sub2ind(size(C),id(k), p(2), cl(2));
                        C(id(k),p(2),cl(2)) = -abs(goodcell);

                        C(id(k),p(1),cl(1)) = goodcell;
                        if sum( abs(C(id(k),p(1),:)) == goodcell) > 1
                            fprintf('Created Duplicate of %d on same probe!! E%dP%d\n',goodcell,e,p(1));
                            errs(end+1) = 6;
                        end
                        dup(id(k),badcell) = -goodcell;
                        if p(1) == p(2)
                            fprintf('Found Duplicate on same probe!! E%dP%d\n',e,p(1));
                            errs(end+1) = 6;
                        end
                        replacelist(end+1,:) = [rid badcell goodcell sum(errs) 0];
                        errlist{size(replacelist,1)} = errs;
                    end
                end
%Check that left over cells with the duplicate are not expts
%where the primary was missing
                if goodcell ~= badcell
                    bid = find(C == badcell);
                    [e,p,cl] = ind2sub(size(C),bid);

                    gid = find(C == goodcell);
                    [ge,gp,gcl] = ind2sub(size(C),gid);
                    
                    [ebid, epid, ecid] = find(C == -badcell);
                    if ~isempty(bid) && ~isempty(gid)
                        comp = CompareCells(DATA, C, Clusters, [goodcell badcell]);
                        if comp.merge > 0 
                            fprintf('Merging %d->%d (merge %d)\n',goodcell, badcell,comp.merge);
                            if ~isempty(comp.nonmatch)
                                cprintf('red','BUT nonmatches at %d\n',sprintf('%d ',comp.nonmatch));
                                errs(end+1) = 7;
                            end
                            [oex, eoid, goid] = intersect(e,ge);
                            if ~isempty(oex)
%This can happen if an earlier replacement (my merging) means that badcell appears in an experiment
% it did not previously seen, so xcorr wast not listed in Groups G at the
% start
                                fprintf('Cells %d and %d both present in E%d\n',oex);
                                for j = 1:length(oex)
                                    efficacy = clust.Efficacy(Clusters,[ge(goid(j)) gp(goid(j)) gcl(goid(j));e(eoid(j)) p(eoid(j)) cl(eoid(j))]);
                                    eff(j) = max(efficacy(1:2));
                                    if eff(j) < 0.2
                                        fprintf('Cell %d in Ex%dP%d.%d not a duplicate of %d\n',badcell,oex(j),p(eoid(j)),cl(eoid(j)),goodcell);
                                        errs(end+1) = 8;
                                    else
                                        errs(end+1) = 9;
                                    end
                                end
                            end
                            C(bid) = goodcell;
                            [aa,bb,cc] = ind2sub(size(C),bid);
                            CC = PC.GetClusterInfo(Clusters,cat(2,aa(:),bb(:),cc(:)),'clst','array');
                            for j = 1:length(bid)
                                replacelist(end+1,:) = [bid(j) badcell goodcell sum(errs) sum(CC{j}.clst == cc(j)+1)];
                                if sum(abs(C(e(j),p(j),:)) == goodcell) > 1
                                    fprintf('Merge Created Duplicate of %d on same probe!! E%dP%d\n',goodcell,e(j),p(j));
                                    errs(end+1) = 9;
                                end
                                errlist{size(replacelist,1)} = errs;
                            end
                            did = find(C == -badcell);
                            C(did) = -goodcell;
                            findcell.SpecialCheck(C, Clusters);
                        else
                            C = findcell.CheckOrphanMerge(DATA, C, Clusters, goodcell, badcell);
                        end
                    end
                    if sum(C(:) == -badcell) && ~sum(C(:) == badcell)
                        oid = find(C == -badcell);
                        fprintf('Setting %d orphaned %d to %d\n',length(oid),-badcell, -goodcell);
                        C(oid) = -goodcell;
                    end
                    if sum(C(:) == goodcell) == 0
                        fprintf('E%dCell%d disappeared\n',goodcell);
                    end
                end
                
            end
            cdup = celllist.FindDuplicates(CellList);
            if length(cdup) > length(lastdup)
                fprintf('New Self Duplicates when seting %d\n',goodcell);
            end

           end
           if nc > 0 %did find some non-self duplicates. Rebuild duplicate list
            [hcount, dup] = CountDuplicates(C ,G, Clusters);
            [a,b] = find(hcount > 0);
            fprintf('Duplicate Listing Still %d duplicates left\n',length(a));
           else
               fprintf('FixDuplicates: Finished marking duplicates first pass\n');
               if length(a)
                   fprintf('Some Self Duplicates  %s %s\n',sprintf(' %d',a),sprintf(' %d',b));
               end
               a = [];
           end
        end
        
%Now run through and check that all duplicates are gone
      for j = 1:length(G)
          if isfield(G{j},'groups')
              gcells = G{j}.cells;
              plist = cat(1,G{j}.xcorrs.p);
          for k = 1:length(G{j}.groups)
              if length(G{j}.groups{k}) > 1
                  cid = G{j}.groups{k};
                  cells = [];
                  for c = 1:length(cid)
                      x = gcells(cid(c));
                      cells(c) = C(x.eid,x.p,x.cl);
                  end
                  ecid = find(cells > 0);
                  dupid = cid(ecid);
                  ec = unique(cells(ecid));
                  if length(ec) > 1
                      fprintf('still duplicates in E%d:',G{j}.expt)
                      a = find(sum(ismember(plist,dupid),2) ==2 & diff(plist,[],2) ~= 0);
                      for ia = 1:length(a)
                          x = G{j}.xcorrs(a(ia));
                          [isdup(ia),score(ia,:)] = PC.IsSameCell(x);
                          fprintf(' (%d) %.1f <-> %.1f %.1f,%.1f', isdup(ia),x.probe(1),x.probe(2),score);
                      end
                      if length(a) > 1
                          fprintf('*');
                      end
                      fprintf('\n');
                      if max(score(:,1) > 0.2)
                      clear nc;
                      for c = 1:length(ec)
                          nc(c) = sum(C(:) == ec(c));
                      end
                      [a,b] = max(nc);
                      goodcell = ec(b);
                      [newlist, Rdetails] = findcell.ResolveDuplicate(DATA, cells, C, Clusters);
                      C = newlist;
                      for c = 1:length(cells)
                          if cells(c) >0 && cells(c) ~= goodcell
                              x = gcells(cid(c));
                              id = sub2ind(size(C),x.eid,x.p,x.cl);
                              oldcell = C(x.eid,x.p,x.cl);
                              [oldp, oldcl] = find(abs(squeeze(C(x.eid,:,:))) == oldcell);
                              oldid = sub2ind(size(C),ones(size(oldp)).*x.eid, oldp, oldcl);
                              C(oldid) = -goodcell;
                              C(x.eid,x.p,x.cl) = -goodcell;
                              findcell.SpecialCheck(C, Clusters);
                              replacelist(end+1,1:5) = [id cells(c) goodcell sum(errs) 0];
                          end
                      end
                      else
                          fprintf('Not relacing becauase < 10%% common spikes\n');
                      end
                  end
              end
          end
          end
      end
      
        [a,b] = find(hcount > 0);
        details.replaced = replaced;
        details.checkpoints = checkpoints; 
        details.replacelist = replacelist;
        details.partiallist = FindPartialDuplicates(C, G, Clusters);
        details.warnings = errlist;


%Now align duplicates so that a cell is on one probe as much as possible
dup = unique(C(find(C < 0)));
for j = 1:length(dup)
    [e,p,cl] = celllist.find(C,-dup(j));
    if isempty(e)
        fprintf('Duplicate %d  has not matching Cells\n',dup(j));
    end
    [a,b] = Counts(p);
    [c,d] = max(a);
    probe = b(d);
    for k = 1:length(e)
        clid = find(squeeze(C(e(k),probe,:)) == dup(j));
        if length(clid) > 1
            fprintf('E%dP%d cell%dDuplicate Cell on One probe\n',e(k),p(k));
            if I(e(k),probe,clid(1)) > I(e(k),probe,clid(2))
                clid = clid(1);
            else
                clid = clid(2);
            end
        end
        if p(k) ~= probe && ~isempty(clid)  %main probe marked as dup
            if (I(e(k),probe,clid) > 4 || I(e(k),probe,clid) > 0.9 .* I(e(k),p(1),cl(1)))  && ...
                    (D(e(k),probe,clid) > 2 || D(e(k),probe,clid) > 0.9 .* D(e(k),p(1),cl(1)))
                C(e(k),probe,clid) = -dup(j);
                findcell.SpecialCheck(C, Clusters);
                C(e(k),p(k),cl(k)) = dup(j);
            end
        end
    end
    cellid = -dup(j);
end

cells = unique(C);
cells = cells(cells>0);
for j = 1:length(cells)
    [e,p] = celllist.find(C,cells(j));
    [a,b] = Counts(e);
    id = find(a > 1);
    for k = 1:length(id)
        eid = b(id(k));
        pid = find(e == eid);
        fprintf('Cell %d Duplicated on E%dP%s\n',cells(j),eid,sprintf('%d,',p(pid)));
        [cp,clid] = find(abs(squeeze(C(eid,:,:))) == cells(j)); %all matches for this expt
        clear vals;
        for c = 1:length(cp)
            A = PC.GetClusterInfo(Clusters,[eid cp(c) clid(c)]);
            vals(c,:) = clust.GetValue(A,0,{'defaultisolation' 'dropi'});
        end
        [ia, ii] = max(vals);
        if ii(1) == ii(2)
            best = ii(1);
        else
            best = ii(1); %best isolation for now. 
        end
        cid = sub2ind(size(C),ones(size(cp)) .* eid, cp, clid);
        C(cid) = -cells(j);
        C(eid,cp(best),clid(best)) = cells(j);
        details.bestduplicate(eid, cells(j)) = cid(best);
    end
end

function SpecialCheck(C)
    if size(C,1) > 43 && C(44,21,1) == C(44,21,2)
        C(44,21,2) = C(44,21,2);
    end

            
function plist =  FindPartialDuplicates(C, G, Clusters)
    plist = [];
    [hcount, dup] = CountDuplicates(C ,G,Clusters,'partial');
    [a,b] = find(hcount > 0);
    for j = 1:length(a)
        [ea, pa, cla] = celllist.find(C, a(j));
        [eb, pb, clb] = celllist.find(C, b(j));
        [aid,ia,ib] = intersect(ea,eb); %expts with both
        dupid = find(dup(:,a(j)) > 0 & dup(:,b(j)) >0); %expts with dup
        for k = 1:length(aid)
            if ismember(aid(k),dupid)
                probe(1) = pa(ia(k)) + cla(ia(k))./10;
                probe(2) = pb(ib(k)) + clb(ib(k))./10;
                xc = PC.FindXCorr(G{aid(k)}.xcorrs,probe);
                [isdup, score] = PC.IsSameCell(xc);
                if isdup == 2
                    plist(end+1,:) = [aid(k) probe];
                end
            end
        end
    end
        
        
function [C, details] = FixDuplicatesV2(DATA, CellList, D, I, Clusters, varargin)
        
        checkpoints = [];
        details.replaced = [];
        cellgroups = {};
        G = getappdata(DATA.toplevel,'AllDuplicates');
        C = CellList;
        [hcount, dup] = CountDuplicates(C ,G, Clusters);
        replaced = [];
        replacelist = [];
        [a,b] = find(hcount > 0);
%Hard to keep track of cells as they get replaced because
%original duplicate listting (G) is based on groupings within
%each expt, not just pairised corrs.
%Could go through G each time a cell is replaced and find all
%matches.  But thats a pain.
       details.duplist = hcount;
        for j = 1:length(a)
            if b(j) == a(j)
            else
                acell = a(j);
                bcell = b(j);
                if ~isempty(replaced)
                    rid = find(replaced(:,1) == a(j));
%beware of incomplete replacesments. replaced(x,y) means SOME
%location of x were set as replacements;
                    if ~isempty(rid)
                        replace(1) = replaced(rid(1),2);
                        acell = replace(1);
                        aleft = sum(C(:) == a(j));
                        if aleft > 0
                            fprintf('Cell%d remains in %d expts despite replacement in %d\n',a(j),aleft,length(rid))
                        end
                    else
                        aleft = 0;
                        replace(1) = 0;
                    end
                    rid = find(replaced(:,1) == b(j));
                    if ~isempty(rid)
                        replace(2) = replaced(rid(1),2);
                        bcell = replace(2);
                        bleft = sum(C(:) == b(j));
                        if bleft > 0
                            fprintf('Cell%d remains in %d expts despite replacement in %d\n',b(j),bleft,length(rid))
                        end
                    else
                        bleft = 0;
                        replace(2) = 0;
                    end
                else
                    replace(1:2) = 0;
                end
                id = find(dup(:,acell) > 0 & dup(:,bcell) >0);
                if isempty(id)
                    if aleft && bleft
                        acell = a(j);
                        bcell = b(j);
                    elseif aleft
                        acell = a(j);
                    elseif bleft
                        bcell = b(j);
                    else
                    end
                    if acell == bcell
                        aleft = sum(C(:) == a(j));
                        bleft = sum(C(:) == b(j));
                        fprintf('%d (%d)->%d(%d) becomes %d->%d\n',a(j),aleft,b(j),bleft,acell,bcell)
                    else                        
                        id = find(dup(:,acell) > 0 & dup(:,bcell) >0);
                    end
                end
                na = sum(C(:) == acell);
                nb = sum(C(:) == bcell);
                clear good;
                clear bad;
                if na < nb
                    good(1:length(id)) = bcell;
                    bad(1:length(id)) = acell;
                    goodcell = bcell;
                    badcell = acell;
                else
                    bad(1:length(id)) = bcell;
                    good(1:length(id)) = acell;
                    goodcell = acell;
                    badcell = bcell;
                end
                    
                fprintf('Cell%d(%d) %d, Cell%d(%d) %d, %d/%d duplicate. Using cell %d Ex%s\n',...
                a(j),replace(1),na,b(j),replace(2),length(id),nb,hcount(a(j),b(j)),goodcell,sprintf(' %d',id));  
                for k = 1:length(id)
                    e = id(k);
                    eC = squeeze(C(id(k),:,:));
                    gid = find(eC == good(k));
                    bid = find(eC == badcell);
                    bid = setdiff(bid,gid); %don't self replace
%A cell might be missing if it was a duplicate with someone
%else and so was reset in a previous loop
                    if ~isempty(bid) && ~isempty(gid)
                        replaced(end+1,1) = badcell;
                        replaced(end,2:5) = [goodcell e a(j) b(j)];
                        if length(gid) > 1
                            fprintf('Duplicate Good cell %d in Expt %d\n',good(k),e)
                            ind2sub(size(eC),gid);
                        end
                        if length(bid) > 1
                            fprintf('Duplicate Bad cell %d in Expt %d\n',bad(k),e)
                            ind2sub(size(eC),gid);
                        end
                        [p(1), cl(1)] = ind2sub(size(eC),gid(1));
                        [p(2:length(bid)+1), cl(2:length(bid)+1)] = ind2sub(size(eC),bid);
                        best(k) = good(k);
                        if D(e,p(2),cl(2)) > D(e,p(1),cl(1)) && I(e,p(2),cl(2)) > I(e,p(1),cl(1))
                            if I(e,p(1),cl(1)) < 4 || D(e,p(1),cl(1)) < 2
                                best(k) = bad(k);
                                bad(k) = good(k);
                                good(k) = best(k);
                                p = fliplr(p);
                                cl = fliplr(cl);
                            end
                        end                       
                        
%anthing currrently recorded as a copy of badcell should be set to a copy of goodcell                        
                        [bpid, bcid] = find(squeeze(C(id(k),:,:)) == -badcell);
                        for c = 1:length(bpid)
                            C(id(k),bpid(c),bcid(c)) = -goodcell;
                        end
%replacelist records original cell# in locattions that get marked as duplicates                        
                        rid = sub2ind(size(C),id(k), p(2), cl(2));
                        replacelist(end+1,:) = [rid badcell goodcell sum(errs)];
                        C(id(k),p(2),cl(2)) = -goodcell;
                        C(id(k),p(1),cl(1)) = goodcell;
                        dup(id(k),badcell) = -goodcell;
                        if p(1) == p(2)
                            fprintf('Found Duplicate on same probe!! E%dP%d\n',e,p(1));
                        end
                    end
                end
%Check that left over cells with the duplicate are not expts
%where the primary was missing
                if goodcell ~= badcell
                    bid = find(C == badcell);
                    [e,p,cl] = ind2sub(size(C),bid);

                    gid = find(C == goodcell);
                    [ge,gp,gcl] = ind2sub(size(C),gid);
                    
                    [ebid, epid, ecid] = find(C == -badcell);
                    if ~isempty(bid) && ~isempty(gid)
                        comp = CompareCells(DATA, C, Clusters, [goodcell badcell]);
                        if comp.merge
                            fprintf('Merging %d->%d\n',goodcell, badcell);
                            [oex, eoid, goid] = intersect(e,ge);
                            if ~isempty(oex)
%This can happen if an earlier replacement (my merging) means that badcell appear in an experiment
% it did not previously seen, so xcorr wast not listed in Groups G at the
% start
                                fprintf('Cells %d and %d both present in E%d\n',oex);
                                for j = 1:length(oex)
                                    efficacy = clust.Efficacy(Clusters,[ge(goid(j)) gp(goid(j)) gcl(goid(j));e(eoid(j)) p(eoid(j)) cl(eoid(j))]);
                                    eff(j) = max(efficacy(1:2));
                                    if eff(j) < 0.2
                                        fprintf('Merge: Cell %d in Ex%dP%d.%d not a duplicate of %d\n',badcell,oex(j),p(eoid(j)),cl(eoid(j)),goodcell);
                                    end
                                end
                            end
                            C(bid) = goodcell;
                            for j = 1:length(bid)
                                replacelist(end+1,:) = [bid(j) badcell goodcell sum(errs)];
                            end
                            did = find(C == -badcell);
                            C(did) = -goodcell;
                        else
                            C = findecll.CheckOrphanMerge(DATA, C, Clusters, goodcell, badcell);
                        end
                    end
                    if sum(C(:) == -badcell) && ~sum(C(:) == badcell)
                        oid = find(C == -badcell);
                        fprintf('Setting %d orphaned %d to %d\n',length(oid),-badcell, -goodcell);
                        C(oid) = -goodcell;
                    end
                end
                
            end
        end
        
%Now run through and check that all duplicates are gone
      for j = 1:length(G)
          if isfield(G{j},'groups')
              gcells = G{j}.cells;
          for k = 1:length(G{j}.groups)
              if length(G{j}.groups{k}) > 1
                  cid = G{j}.groups{k};
                  cells = [];
                  for c = 1:length(cid)
                      x = gcells(cid(c));
                      cells(c) = C(x.eid,x.p,x.cl);
                  end
                  ec = unique(cells(cells>0));
                  if length(ec) > 1
                      fprintf('still duplicates in E%d:%s\n',G{j}.eid,sprintf(' %d',ec));
                      clear nc;
                      for c = 1:length(ec)
                          nc(c) = sum(C(:) == ec(c));
                      end
                      [a,b] = max(nc);
                      goodcell = ec(b);
                      for c = 1:length(cells)
                          if cells(c) >0 && cells(c) ~= goodcell
                              x = gcells(cid(c));
                              id = sub2ind(size(C),x.eid,x.p,x.cl);
                              C(x.eid,x.p,x.cl) = -goodcell;
                              replacelist(end+1,:) = [id cells(c) goodcell sum(errs)];
                          end
                      end
                  end
              end
          end
          end
      end
      
        [a,b] = find(hcount > 0);
        details.replaced = replaced;
        details.checkpoints = checkpoints; 
        details.replacelist = replacelist;
%Now align duplicates so that a cell is on one probe as much as possible
dup = unique(C(find(C < 0)));
for j = 1:length(dup)
    [e,p,cl] = celllist.find(C,-dup(j));
    if isempty(e)
        fprintf('Duplicate %d  has not matching Cells\n',dup(j));
    end
    [a,b] = Counts(p);
    [c,d] = max(a);
    probe = b(d);
    for k = 1:length(e)
        clid = find(squeeze(C(e(k),probe,:)) == dup(j));
        if length(clid) > 1
            fprintf('E%dP%d cell%dDuplicate Cell on One probe\n',e(k),p(k));
            if I(e(k),probe,clid(1)) > I(e(k),probe,clid(2))
                clid = clid(1);
            else
                clid = clid(2);
            end
        end
        if p(k) ~= probe && ~isempty(clid)  %main probe marked as dup
            if (I(e(k),probe,clid) > 4 || I(e(k),probe,clid) > 0.9 .* I(e(k),p(1),cl(1)))  && ...
                    (D(e(k),probe,clid) > 2 || D(e(k),probe,clid) > 0.9 .* D(e(k),p(1),cl(1)))
                C(e(k),probe,clid) = -dup(j);
                C(e(k),p(k),cl(k)) = dup(j);
            end
        end
    end
end

    
   function [C, cells] = ReNumberCellList(C)
    
       cells = unique(C(C>0));
       for j = 1:length(cells)
           id = find(C == cells(j));
           C(id) = j;
           id = find(C == -cells(j));
           C(id) = -j;
       end
            
       
    function CompareSpikeShape(A,B, varargin)
        c = mycolors;
        hold off;
        Sa = A.MeanSpike.ms;
        Sb = B.MeanSpike.ms;
        t = sprintf('E%dP%dcl%d vs E%dP%dvs%d',A.exptno,A.probe,A.cluster,B.exptno,B.probe,B.cluster);
        for j = 1:size(Sa,1)
            plot(Sa(j,:),'-','color',c{j});
            hold on;
            plot(Sb(j,:),'--','color',c{j});
        end
        title(t);


        
function Clusters = FindDoubleTriggers(Clusters)
    for j = 1:length(Clusters)
        for k = 1:length(Clusters{j});
            C = Clusters{j}{k};
            cls = unique(C.clst);
            for cj = 1:length(cls)
                for ck = 1:cj-1
                    ta = C.t(C.clst==cj);
                    tb = C.t(C.clst==ck);
                    [e, details] = CalcEfficacy(ta,tb);
                    nspk = details.nspk;
                    if e(1) > 0.3
                        GetFigure('xcorr');
                        fprintf('Repeat Cluster E%dP%d: %d(%d)-%d(%d) %.2f',j,k,cj,nspk(1),ck,nspk(2),e(1));
                        t = [-0.002:0.0001:0.002];
                        xc = xcorrtimes(ta,tb,'times',t);
                        plot(t(2:end),xc);
                        title(sprintf('cl%d(%d)->cl%d(%d)\n',cj,nspk(1),ck,nspk(2)));
                        if ck > 1
                           if nspk(1) < nspk(2)
                               badcl = cj;
                           else
                               badcl = ck;
                           end
                           fprintf('Cl%d should go',badcl);
                        else
                            
                        end
                        fprintf('\n');
                    end
                end
            end           
        end
    end
            
    
function swaps = ListClusterSwaps(C)

        swaps = {};
        for j = 1:length(C)
            for k = 1:length(C{j})
                CC = C{j}{k};
                if isfield(CC,'swapped')
                    swaps{j,k} = CC.swapped;
                end
            end
        end

        
function PlotCellQ(DATA, cellq)
        


GetFigure('CellLines','parent',DATA.toplevel);
colors = mycolors(10);
symbols = 'osx^osx^osx^osx^osx^osx^';
hold off;
for j = 1:length(cellq)
    X = cellq(j);
    z = rem(j,10)./50; %small offset for visibility
    pos = cat(1, [X.expt X.probe X.cl], X.matches.pos);
    [~,id] = sort(pos(:,1));
    pos = pos(id,:);
    c = colors{1 + mod(j-1,10)};
    h(j) = plot(z+pos(:,2)+pos(:,3)./10,pos(:,1),'-','marker',symbols(ceil(j/10)),'color',c,'buttondownfcn',{@PC.HitImage, 'lineplot' 'autolist'});
    linedata.cell = j;
    linedata.start = [X.expt X.probe X.cl];
    set(h(j),'UserData',linedata);
    hold on;
    text(X.probe+(X.cl-1)/10,X.expt,sprintf('%d',j),'color',c,'horizontalalignment','right');
end
set(gca,'ydir','reverse');
    
function CheckDetails(FX, varargin)
    
checks = {};
for j = 1:length(varargin)
    if ischar(varargin{j})
        checks{end+1} = varargin{j};
    end
    j = j+1;
end

if isempty(checks)
    checks = 'alignfit';
end
    
for c = 1:length(FX.cellq)
    M = FX.cellq(c).matches;
    for j = 1:length(M)
        if ~isempty(M(j).aligndata)
            A = M(j).aligndata;
            type = M(j).aligndata.type;
            if strcmp(type,'alignfit') && cellstrcmp(type,checks);
                fprintf('E%dP%d.%d %s',M(j).pos(1:3),type);
                fprintf(' E%dP%d.%d Fix%d',A.pos(1,:),A.usefix);
                fprintf('\n');
            end
            if strcmp(type,'mahalmatrix') && cellstrcmp(type,checks);
                fprintf('E%dP%d.%d %s',M(j).pos(1:3),type);
                fprintf(' P%d.%d',A.probes(1,:));
                fprintf('\n');
            end
        elseif cellstrcmp('listall',varargin)
            fprintf('E%dP%d.%d\n',M.pos(1:3));
        end
    end
end

function ApplyFixes(DATA)

    ts = now;
    FX = getappdata(DATA.toplevel,'AutoCellList')
    Clusters = PC.GetValue(DATA.toplevel,'Clusters');
    for j = 1:length(FX.fixlist)
        cpos = FX.fixlist(j,:);
        eid = find(DATA.exptlist == cpos(1));
        if length(eid) == 1
            fprintf('Fixing E%dP%d  with Fit %.1f. Requested by %d\n',cpos);
            res = PC.CallAllVPcs(DATA.toplevel, eid,cpos(2),'setfit',cpos(3),'saveclusters');
        end
    end
    mytoc(ts);

function mem = MemoryUse(DATA)
    
    figs = findobj(allchild(0),'flat','type','figure');
    for j = 1:length(figs)
        f = floor(double(figs{j}));
        app{f} = gui.FigMemoryUse(figs(j));
        mem(f).tag = app{f}.tag;
        mem(f).total = app{f}.total;
    end

    

    
        
function LogState(DATA, finddetails)
    try
    if isfield(finddetails,'stage') && finddetails.stage > 0
        state.stage = finddetails.stage;
        state.ncells = length(finddetails.cellq);
    else
        state.stage = -1;
        state.ncells = 0;
    end
    [a,b] = memory;
    gb=(1024.*1024.*1024);
    state.mem(1) = b.SystemMemory.Available./gb;
    state.mem(2) = b.PhysicalMemory.Available./gb;
    state.mem(3) = a.MemUsedMATLAB./gb;
    mem = MemoryUse(DATA);
    state.figmemory = [mem.total];
    state.mem(4) = sum(state.figmemory);
    state.date = now;
    if isfield(DATA,'autodatadir')
        logname = [DATA.autodatadir '/FindLog.mat'];
    end
    X = my.load(logname,'safe');
    if isfield(X,'states')
        states = X.states;
        states(end+1) = state;
    else
        states(1) = state;
    end
    save(logname,'states');
    end
    