function [FX, allcheck] = CheckAlignment(DATA, FX, varargin)
%PC.CheckAlignment(DATA, finddata,...) Checks cluster shape is consistent
%across expt boundaries
%PC.CheckAlignment(DATA) uses current appdata for the AutoCellList
%PC.CheckAlignment(DATA,'list')
%               ....,'cells',cellids) only do checks on cells in list
cellids = [];
applyfixes = 0;
boundaryonly = 0;
DATA = GetDataFromFig(DATA);
if nargin < 2 || ischar(FX)
    if nargin > 1
        varargin = {FX varargin{:}};
    end
    FX = getappdata(DATA.toplevel,'AutoCellList');
end

usenewclusters = 0;
checkgaps = 1;
%ellipse fitting can produce Clusters with bad clsts 
%and is time consuming.  Not for now....
args = {'noreverse' 'noellipse'};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'apply',5)
        applyfixes = 1;
        if strcmpi(varargin{j},'applytest')
            applyfixes = 2;
        end
    elseif strncmpi(varargin{j},'boundary',5)
        boundaryonly = 1;
    elseif strncmpi(varargin{j},'cells',5)
        j = j+1;
        cellids = varargin{j};
    elseif strncmpi(varargin{j},'list',4)
        PrintAlignErrs(FX, varargin{:});
        return;
    elseif strncmpi(varargin{j},'newclusters',8)
        usenewclusters = 1;
    elseif strncmpi(varargin{j},'nogaps',5)
        checkgaps = 0;
    elseif strncmpi(varargin{j},'reapply',5)
        ApplyFixes(DATA, FX);
        return;
    elseif strncmpi(varargin{j},'watch',5)
        FX.state.watchprogress = 1;
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
    
end


changes = {};
fixes = {};
cellchanges = [];
if isfield(DATA.autolist,'CellList')
    CellList = DATA.autolist.CellList;
elseif ~isfield(FX,'CellList')
    CellList = FX.CellLists{end};
else
    CellList = FX.CellList;
end
if applyfixes == 1
    setappdata(DATA.toplevel,'AutoCellListBackup',CellList);
end
 cells = unique(CellList);
 cells = cells(cells>0);
 if ~isempty(cellids)
     cells = intersect(cells, cellids);
 end
 if isappdata(DATA.toplevel,'NewClusters') && usenewclusters
     Clusters = getappdata(DATA.toplevel,'NewClusters');
 else
     Clusters = PC.GetValue(DATA,'Clusters','withfits');
 end
 allcheck = [];
FX.state.toplevel = DATA.toplevel;
FX.state = AddField(FX.state,{'verbose','watchprogress'},1);
 setappdata(DATA.toplevel,'PreClusters',Clusters);
 %Fill gaps before checking all transitions
 [~, errs] = celllist.Check(CellList,Clusters);
 if isfield(errs,'type') && sum(strcmp('clusters',{errs.type}))
     id = find(strcmp('clusters',{errs.type}));
 end
 if checkgaps
     [gaptransitions, newCellList, Clusters] = CheckBoundaryLosses(cells, CellList, Clusters, FX,applyfixes);
     findcell.SetDetails(DATA, FX, newCellList, FX.stage, FX.stage);
     FX.gaptransitions = CellToStruct(gaptransitions);
     if isempty(cellids) %only setappdata if all cells were done
         setappdata(DATA.toplevel,'AutoCellList',FX);
     end
 end
 if isfield(errs,'type') && sum(strcmp('clusters',{errs.type}))
     id = find(strcmp('clusters',{errs.type}));
 end
 if boundaryonly
     badtransitions = CheckBoundaryMahal(cells, CellList, Clusters, FX);
     if ~isempty(badtransitions)
         FX.badtransitions = CellToStruct(badtransitions);
         if isempty(cellids)
             setappdata(DATA.toplevel,'AutoCellList',FX);
         end
     end
     return;
 end
 ts = now;
 for j = 1:length(cells)
     set(FX.state.toplevel,'Name',sprintf('Cell %d Checking Alignemnt',cells(j)));
     drawnow;
     [eid, pid, cid] = celllist.find(CellList, cells(j));
     ediff = diff(eid)-1;
     gaps = find(ediff > 0 & ediff < 3);
     missing = [];
     for g = 1:length(gaps)
         for k = 1:ediff(gaps(g))
             missing(end+1,:) = eid(gaps(g))+k;
             [~,id] = min(abs(missing(end)-eid));
             eid(end+1) = missing(end);
             pid(length(eid)) = pid(id);
             cid(length(eid)) = cid(id);
         end
     end
     if ~isempty(missing)
         [eid,id] = sort(eid);
         pid = pid(id);
         cid = cid(id);
     end
     fprintf('\n\nChecking Cell %d\n',cells(j));
     for e = 1:length(eid)-1
         pos = [eid(e) pid(e) 1; eid(e+1) pid(e+1) 1];
         na = sum(CellList(eid(e),pid(e),:) == cells(j));
         nb = sum(CellList(eid(e+1),pid(e+1),:) == cells(j));
         C = PC.GetClusterInfo(Clusters,pos,'clst');
         [a,b] = clust.AlignAutoFits(C,'checkcl',cid(e),args{:});
         b.cell = cells(j);
         allcheck(end+1).cell = b.cell;
         allcheck(end).aligned = b.aligned;
         b.stage = 6;
         fixed = b;
         fixed.good = 0;
         if ismember(b.aligned,[0 3])
             b.flipped = 0;
             [changes, fixes, b,Clusters, CellList] = Align(DATA,Clusters, changes,fixes, b,applyfixes, CellList);
             b.pos = pos;
             fixed.good = 1;
         elseif b.aligned == -2  %cid(e) isn't in Cluster
             fixed.fail = 'missingcluster';
             fprintf('%s Does not have Cluster %d\n',str(C{1}),cid(e));
         elseif b.aligned == -1  %cid(e) isn't in Cluster
             fixed.fail = 'missingcluster';
             fprintf('%s Does not have Cluster %d\n',str(C{2}),cid(e+1));
         elseif b.match ~= cid(e+1) %cell is not best match
             [M,D] = clust.MahalMatrix(C,'plotmatch',cid(e));
             allcheck(end).match = [b.match, D.bestmatch];
             if D.bestmatch == b.match %agree its wrong
                 if CellList(eid(e+1),pid(e+1),b.match) <= 0
                     CellList(eid(e+1),pid(e+1),b.match) = cells(j);
                     CellList(eid(e+1),pid(e+1),cid(e+1)) = 0;
                 else
                     ocell= CellList(eid(e+1),pid(e+1),b.match);
                     [a,c] = clust.AlignAutoFits(C([2 1]),'checkcl',cid(e+1),args{:});
                     [bM,bD] = clust.MahalMatrix(C([2 1]),'plotmatch',cid(e+1));
                     if c.match == bD.bestmatch
                         if CellList(eid(e),pid(e),c.match) == ocell %simple swap
                             fprintf('Swapping Cells %d <->%d in %s\n',cells(j),ocell,str(C{2}));
                             oldcell = CellList(eid(e),pid(e),c.match);
                             CellList(eid(e+1),pid(e+1),b.match) = cells(j);
                             CellList(eid(e+1),pid(e+1),cid(e+1)) = oldcell;
                             cid(e+1) = b.match; %so that next loop used new cell
                         else
                             fixed.fail = 'swap';
                             fprintf('Cant resolve cell swapping at %s,%s\n',str(C{1}),str(C{2}));
                         end
                     else
                             fixed.fail = 'matchtwo';
                             fprintf('XC and mahl differetn matches at %s,%s\n',str(C{2}),str(C{1}));                         
                     end
                 end
             else
                 fixed.fail = 'match';
                 fprintf('XC and mahl differetn matches at %s,%s\n',str(C{1}),str(C{2}));                 
             end
         end
         if fixed.good == 0
             changes{end+1} = rmfield(fixed,'good');
             fixes{end+1} = fixed;
         end
         if isfield(b,'matches') && sum(ismember(b.matches,cid(e))) == 0
             [a,b] = clust.AlignAutoFits(C([2 1]),'checkcl',b.matches(1),args{:});
         end
         [a,b] = clust.AlignAutoFits(C([2 1]),'checkcl',cid(e+1),args{:});
         b.cell = cells(j);
         b.stage = 7;
         if ismember(b.aligned,[0 3])
             b.pos = flipud(pos);
             b.flipped = 1;
             [changes, fixes, b,Clusters, CellList] = Align(DATA, Clusters,changes,fixes, b,applyfixes, CellList);
         end
     end
     if isfield(errs,'type') && sum(strcmp('clusters',{errs.type}))
         id = find(strcmp('clusters',{errs.type}));
     end
 end 
 set(DATA.toplevel,'Name',DATA.datadir);

 for j = 1:length(changes)
     if isfield(changes{j},'fixlist') 
         f = changes{j}.fixlist;
         if ~isempty(f)
             FX.fixlist(end+1,1:length(f)) = f;
             FX.fixfit(f(1),f(2)) = f(3);
         end
     end
 end

X.duration = mytoc(ts);
fprintf('Took %.1f sec\n',X.duration);
if applyfixes == 1
    listtype = PC.GetValue(DATA, 'listtype');
    if strcmp(listtype,'autolist')
        setappdata(DATA.toplevel,'AutoClusters',Clusters);
    else
        setappdata(DATA.toplevel,'Clusters',Clusters);
    end
    outname = [DATA.autodatadir '/FixedClusters.mat'];
    save(outname,'Clusters');
elseif applyfixes ==2
    setappdata(DATA.toplevel,'NewClusters',Clusters);    
end

badtransitions = CheckBoundaryMahal(cells, CellList, Clusters, FX);
if ~isempty(badtransitions)
    FX.badtransitions = CellToStruct(badtransitions);
end
 if exist('fixes','var')
     X.fixes = CellToStruct(fixes);
     X.changes = CellToStruct(changes);
     FX.alignerrs = X.changes;
     FX.alignfixes = X.fixes;
     setappdata(DATA.toplevel,'AutoCellList',FX);
 end
 

   
function [usematch, details] = FindMatch(cl, M, varargin)
%given matrix of mahal distances, find best match for cell cl 
   if isempty(cl)
       usematch = 0;
       return;
   end
   d = reshape(M(cl,:,:),size(M,2),size(M,3));
   [~,matches] = min(d,[],1);
   [~,matches(end+1)] = min(min(d'));
   [~,matches(end+1)] = min(max(d'));
   [~,matches(end+1)] = min(mean(d'));
   [a,b] = Counts(matches);
   usematch = 0;
   [maxcount,id] = max(a);
   if length(a) ==1
       usematch = b;
   elseif max(a) > 4 %5/7 matches the same
       usematch = b(id);
   else
       details.badreason = 'nomatch';
   end
   details.bestmatch = b(id);
   if usematch       
       %found a match. Check for uniqueness
       details.matchd = squeeze(M(cl,usematch,:));
       if ~cellstrcmp('recurse',varargin);
           for j = 1:size(M,1)
               [allmatches(j), X{j}] = FindMatch(j,M,'recurse');
               if isfield(X{j},'matchd')
                   details.allmatchd(j) = prctile(X{j}.matchd,40);
               end
           end
           details.allmatches = allmatches;
           xmatch = setdiff(find(allmatches>0),cl);
           if ~isempty(xmatch) %more than one match
               nmatchi = find(allmatches == usematch);
               nmatch = length(nmatchi);
               X = CellToStruct(X);
               if sum(details.allmatchd(nmatchi) < 4.5) > 1
                   usematch = 0;
                   details.badreason = 'twomatches';
               end
           end
          
       end
       if sum(details.matchd < 4.5) < 2
           usematch = 0;
           details.badreaon = 'distance';
       end
   end
   
        
function [gaptransitions, CellList, Clusters] = CheckBoundaryLosses(cells, CellList, Clusters, FX, applyfixes)
 gaptransitions = {};
DATA = get(FX.state.toplevel,'UserData');        
errs = {};
fixpos = [];
for j = 1:length(cells)
     [eid, pid, cid] = celllist.find(CellList, cells(j));
     [xeid, xpid, xcid] = celllist.find(CellList, -cells(j));
     [xeid,xid] = setdiff(xeid,eid);
     xeid = unique(xeid);
     for e = 1:length(xeid)
         [CellList, details{e}] = findcell.BestDuplicate(CellList, Clusters, xeid(e),cells(j));
         if isfield(details{e},'pos')
             fixpos(end+1,:) = [details{e}.pos cells(j)];
         end
     end     
end

if ~isempty(fixpos)
    errs{1}.fixpos = fixpos;
end

for j = 1:length(cells)
     set(FX.state.toplevel,'Name',sprintf('Cell %d Checking Gaps',cells(j)));
     drawnow;
     [eid, pid, cid] = celllist.find(CellList, cells(j));
         
     ediff = diff(eid)-1;

     gaps = find(ediff > 3);
%For larger gaps, need to follow the cluster     
     for e = 1:length(gaps)
         expts = eid(gaps(e)):eid(gaps(e)+1)-1;
         if sum(ismember(expts,FX.shakelist)) == 0
         for k = 1:length(expts)-1
             pos = [expts(k) pid(e) 1; expts(k)+1 pid(e) 1];
             [CellList, fix, Clusters] = findcell.CheckGap(DATA, FX, pos, Clusters, CellList,cells(j),applyfixes);
             errs{end+1} = fix;
             if ~isfield(fix,'fixlist')
             end
         end
         end
     end
     maxe = min([FX.shakelist size(CellList,1)]);
     e = eid(end);
     while e < maxe
         pos = [e pid(end) 1; e+1 pid(end) 1];
         [newCellList, fix, newClusters] = CheckGap(DATA, FX, pos, Clusters, CellList,cells(j),applyfixes);
         Clusters = newClusters;
         CellList = newCellList;
         errs{end+1} = fix;
         if isfield(fix,'fixlist') || fix.usematch > 0
             e = e+1;
         else
             e = maxe;
         end
     end
     
     gaps = find(ediff > 0 & ediff < 3);
     missing = [];
     nx =0;
     allpos = [];
     for gi = 1:length(gaps)
         e = gaps(gi);
         pos = [eid(e) pid(e) cid(e); eid(e)+1 pid(e) 1];
         if isempty(allpos)
             allpos(:,:,1) = pos;
         else
             allpos(:,:,end+1) = pos;
         end
         if ediff(e)  == 2
             pos = [eid(e+1) pid(e+1) cid(e+1); eid(e+1)-1 pid(e+1) 1];
             allpos(:,:,end+1) = pos;
         end
     end
     if isempty(allpos)
         nx = 0;
     else
         nx = size(allpos,3);
     end
     for e = 1:nx
         pos = allpos(:,:,e);
         [newCellList, fix, newClusters] = findcell.CheckGap(DATA, FX, pos, Clusters, CellList,cells(j),applyfixes);
         CellList = newCellList;
         Clusters = newClusters;
         if isfield(fix, 'pos')
             errs{end+1} = fix;
         end
     end
 end
 gaptransitions = errs;
set(FX.state.toplevel,'Name',DATA.datadir);

function badtransitions = CheckBoundaryMahal(cells, CellList, Clusters, FX)
        %Check Boundary transitions after making reclusters above 
 
DATA = get(FX.state.toplevel,'UserData');        
 badtransitions = {};
 for j = 1:length(cells)
     [eid, pid, cid] = celllist.find(CellList, cells(j));
     ediff = diff(eid)-1;
     gaps = find(ediff > 0 & ediff < 3);
     missing = [];
     for e = 1:length(eid)-1
         pos = [eid(e) pid(e) 1; eid(e+1) pid(e+1) 1];
         aid = find(CellList(eid(e),pid(e),:) == cells(j));
         aid = aid(1);
         bid = find(CellList(eid(e+1),pid(e+1),:) == cells(j));
         bid = bid(1); %in case cell is there twice
         C = PC.GetClusterInfo(Clusters,pos,'clst');
         [M, D] = clust.MahalMatrix(C);
         [~, badida] = CheckClusterMatch(Clusters, CellList, pos(1,:));
         ok(1) = ~sum(badida>0)>0;
         [~, badid] = CheckClusterMatch(Clusters, CellList, pos(2,:));
         ok(2) = ~sum(badid>0)>0;
         D.clusterok = ok;
         D.M = M;
         D.pos = pos;
         D.cell = cells(j);
         D.shake =  sum(ismember(eid(e):eid(e+1),FX.shakelist));
         if sum(ok) < 2
             D.err = unique([badid badida]);
             badtransitions{end+1} = D;
         elseif ~isempty(aid) && ~isempty(bid)
             a = find(D.cls{1} == aid(1)+1);
             b = find(D.cls{2} == bid(1)+1);
             md = squeeze(M(a,b,:)); %all 4 mahal distnaces betweeen these two clusters
             c = minmax(md);
             if md(2) < 4.5 && md(3) < 4.5
                 goodmatch = 3;
             elseif md(1) < 4.5 && md(4) < 4.5
                 goodmatch = 2;
             elseif sum(md <4.5) > 1
                 goodmatch = 5;
             else 
                 goodmatch = 0;
             end
             minM = squeeze(min(M,[],3));
             maxM = squeeze(max(M,[],3));
%ai and bi are best macthes For C1 in C2             
             [~, ai] = min(minM(aid,:));
             [~, bi] = min(maxM(aid,:));
             if sum([ai bi] ~= bid)
                 if (ai == bi) % best is same
                     bettermatch = ai;
                     othermatch = bi;
                 elseif bi == aid;
                     bettermatch = ai;
                     othermatch = bi;
                 else
                     bettermatch = bi;
                     othermatch = ai;
                 end
             else
                 bettermatch = 0;
                 othermatch = 0;
             end
             if FX.state.verbose
                 fprintf('E%dP%d.%d - E%dP%d.%d',eid(e),pid(e),aid,eid(e+1),pid(e+1),bid);
                 fprintf('%7.2f',squeeze(M(a,b,:)));
                 fprintf('Shake %d Cell%d',D.shake,cells(j));
                 fprintf('\n');  
             end
             if FX.state.watchprogress
                 cmp = PC.CompareQuickSpks(FX.state.toplevel, pos);
                 figure(cmp.fig(1));
                 text(2,1,sprintf('Cell%d',cells(j)),'fontsize',18,'color',DATA.colors{aid+1});
                 figure(cmp.fig(2));
                 text(2,1,sprintf('Cell%d',cells(j)),'fontsize',18,'color',DATA.colors{bid+1});
                 drawnow;
             end
             if max(md) > 20 || goodmatch == 0 || bettermatch > 0
                 [xc, xcdetail] = PC.ShapeCorr(C{1},C{2},'delays',5,'proberange',FX.corrproberange,'drifts',1,FX.ArrayConfig);
                 D.err = 0;
                 fprintf('E%dP%d.%d - E%dP%d.%d',eid(e),pid(e),aid,eid(e+1),pid(e+1),bid);
                 fprintf('%7.2f',squeeze(M(a,b,:)));
                 fprintf(' probeshift %d Shake %d',xcdetail.probeshift,D.shake);
                 fprintf('\n');
                 D.M = M;
                 D.pos = pos;
                 D.mahal = md;
                 D.cid = [a b];
                 D.cell = cells(j);
                 D.oldmatch = bid(1);
                 D.goodmatch = goodmatch;
                 D.bettermatch = [bettermatch othermatch];
                 D.probeshift = xcdetail.probeshift;
                 badtransitions{end+1} = D;
               
             end
         end
     end
  end
 
function ok = GoodFix(fix, CellList)

ok = 1;
%maybe fix.matchdistance should only check cells rather than max?
if ~isfield(fix,'matchdistance') 
    ok =0;
elseif max(fix.matchdistance) > 4
    ok = 0;
    if min(fix.matchdistance) < 4
%fix.duplicates are the clusters in C{2} that match C{1}
%if we are splitting C{1}, then only matters that the match is good
%for cells
%fix.
        e = fix.pos(2,1);
        p = fix.pos(2,2);
        cells = squeeze(CellList(e,p,fix.duplicates)); %cells in C{2}
        cls = fix.duplicates(fix.matchdistance > 4); %bad matches
        if max(fix.matchdistance(cells ~= 0)) < 4
            ok = 1;
        end
    end
elseif ~isfield(fix,'usefix') 
    ok = 0;
elseif fix.usefix == 0
    ok = 0;
elseif fix.usefix == 3 %not ready to use these yet
    ok = 0;
elseif fix.usefix == 1 %special checks for collapsing C{2}
end


function [changes, fixes, b, Clusters, CellList] = Align(DATA, Clusters, changes, fixes, b, applyfixes, CellList)
   changes{end+1} = b;
   if b.bestfix || b.usefix
       fprintf('Fix(%d) %d %d\n',b.flipped,length(changes),b.usefix);
       if b.usefix
           fixes{length(changes)} = b.newfit(b.usefix); %should have fitnumber
           fixes{length(changes)}.fitused = b.usefix;
       else
           fixes{length(changes)} = b.newfit(b.bestfix);
           fixes{length(changes)}.fitused = b.bestfix;
       end
   else
       fprintf('Fix(%d) %d Not usuable\n',b.flipped,length(changes));
       fixes{length(changes)}.fitnumber = 0;
       fixes{length(changes)}.fitused = 0;
   end
   if applyfixes
       dofix = clust.GoodFix(changes{end}, CellList);
   else
       dofix = 0;
   end
   if dofix > 0
       changes{end}.fixed = changes{end}.usefix;
       [Clusters, CellList, D] = findcell.ApplyFix(DATA, Clusters, fixes{end}, changes{end}, CellList);
       changes{end} = CopyFields(changes{end},D,{'fixlist' 'fixerr'});

       fixed(D.e,D.p) = changes{end}.usefix;
       if isfield(D,'changecell')
           cellchanges(length(changes),:) = D.changecell;
       end
       if isfield(D,'celltype')
           listfix{D.e,D.p} = D.celltype;
       end
       if isfield(D,'warning')
           changes{end}.warning = D.warning;
           changes{end}.warncode = 1;
           cprintf('red','%s\n',D.warning);
       end
   else
       changes{end}.fixed = dofix;
   end

   
function [Clusters] = ReApplyFix(DATA, Clusters, fix)
%Apply a fix from a fixlist
%fix.pos has only the position of the cluster needing fixing.
%don't need to work on CellList. This should be making Clusters Match That.
   D.e = find(DATA.exptlist == fix(1));
   D.p = fix(2);
   C = PC.GetClusterInfo(DATA, [D.e D.p 1]);

   res = PC.CallAllVPcs(DATA, D.e,D.p,'setfit',fix(3),'noninteractive');
   cpfields = clust.autofitfields;
   Clusters = clust.Copy(Clusters, res.cluster, [D.e D.p],cpfields);


function [ok,  badcells] = CheckClusterMatch(Clusters, CellList, pos)

badcells = [];
C = PC.GetClusterInfo(Clusters,pos);
aid = find(CellList(pos(1),pos(2),:) ~= 0);
cls = unique(C.clst)-1;
ok = 1;
okid = ismember(aid,cls);
if sum(~ismember(aid,cls))
    ok = 0;
    id = setdiff(aid,cls);
    fprintf('%s Doesnt have Cl%d for Cell %d\n',str(C),id(1),CellList(pos(1,1),pos(1,2),id(1)));
    badcells(id) = CellList(pos(1,1),pos(1,2),id);
end

function ApplyFixes(DATA, FX, varargin)

Clusters = PC.GetValue(DATA.toplevel,'Clusters');
CellList = PC.GetValue(DATA.toplevel,'CellList');

if isfield(FX,'fixlist')
    fixlist = FX.fixlist;

    if isfield(FX,'gaptransitions')
        gapfix = [];
        for j = 1:length(FX.gaptransitions)
            if ~isnan(FX.gaptransitions(j).fixlist(1))
                gapfix(end+1,:) = FX.gaptransitions(j).fixlist;
            end
        end
        fixlist = cat(1,gapfix, fixlist);
    end
    exs = fixlist(:,1) + fixlist(:,2)./100;
    [uex, eid] = unique(exs,'last');
    for j = 1:length(eid)
        Clusters= ReApplyFix(DATA, Clusters, fixlist(eid(j),:), varargin{:});
    end
else
    for j = 1:length(FX.alignerrs)
        dofix(j) = clust.GoodFix(FX.alignerrs(j), CellList);
        if dofix(j) > 0
            [Clusters, CellList, D] = findcell.ApplyFix(DATA, Clusters, FX.alignfixes(j), FX.alignerrs(j), CellList, varargin)
        end
    end
end
setappdata(DATA.toplevel,'Clusters',Clusters);
setappdata(DATA.toplevel,'AutoClusters',Clusters);
outname = [DATA.datadir '/FixedClusters.mat'];
save(outname,'Clusters','-v7.3');
%Applies fixes to clusters, then updates clsuters strcuture. 
     
function PrintAlignErrs(FX, varargin)

    expts = 1:100;
    probes = 1:100;
    listfix = 0;

j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'expts',4)
        j = j+1;
        expts = varargin{j};
    elseif strncmp(varargin{j},'listfix',7)
        listfix = 1;
    elseif strncmp(varargin{j},'probes',4)
        j = j+1;
        probes = varargin{j};
    end
    j = j+1;
end

for j = 1:length(FX.alignerrs)
    A = FX.alignerrs(j);
    goode = sum(ismember(A.pos(:,1),expts));
    goodp = sum(ismember(A.pos(:,2),probes));
    if goode && goodp
        showpoint = 1;
        if listfix 
            if A.usefix == 0 && listfix ==1
                showpoint = 0;
            end
            if max(A.matchdistance) > 7
                showpoint = 0;
            end
        end
    else
        showpoint = 0;
    end
    if showpoint
    if size(A.pos,1) > 1
        fprintf('%d %d.%d <-> %d %d.%d: ',A.pos(1,1:3),A.pos(2,1:3));
    else
        fprintf('%d %d.%d',A.pos(1,1:3));
    end
    fprintf('Fix%d',A.usefix);
    if isfield(A,'matchdistance') && length(A.matchdistance) >1
        fprintf('  Mahal %.2f %.2f',A.matchdistance(1:2));
    end
    fprintf(' #%d',j);
    fprintf('\n');
    end
end