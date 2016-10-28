function [fit, details] = CombineFits(DATA, fits,varargin)
%AllV.CombineFits(DATA, fits)  match clusters across
%fits, and decide on best final grouping
plottype = 'none';
slimfits = 0;
slimfit = 0;
consensus = 0;
doall = 0;

j = 1;
while j <= length(varargin)
    if strcmpi(varargin{j},'all')
        doall = 1;
    elseif strcmpi(varargin{j},'expandall')
        doall = 2;
    elseif strcmpi(varargin{j},'consensus')
        consensus = 1;
    elseif strcmpi(varargin{j},'image')
        plottype = 'image';
    elseif strcmpi(varargin{j},'slim')
        slimfits = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            slimfit = varargin{j};
        end
    end
    j = j+1;
end

for j = 1:length(fits)
    if ~isfield(fits{j},'usestdtemplate')
        if j == 8
            fits{j}.usestdtemplates = 0;
        else
            fits{j}.usestdtemplates = 1;
        end
    end
end

if doall == 2
    nf = 0;
    for j = 1:length(fits)
        nf = nf+1;
        sfits{nf} = fits{j};
        details.fiti(nf,:) = [j 0];
        if isfield(fits{j},'alternates')
            for k = 1:length(fits{j}.alternates)
                nf = nf+1;
                sfits{nf} = AllV.SetFit(fits{j},fits{j}.alternates{k}.type);
                details.fiti(nf,:) = [j k];
            end
        end
    end
    details.groups = GroupFits(sfits);
    for j = 1:length(details.groups)
        for k = 1:length(details.groups(j).group)
            id = details.groups(j).group(k);
            sfits{id}.group = j;
        end
    end
    fit = sfits;
    return;
elseif doall == 1
    if DATA.usestdtemplates == 0
        DATA.usestdtemplates =1;
        DATA = AllV.CalcTemplateScores(DATA);
    end
    for j = 1:length(fits)
        if isfield(fits{j},'Template')
            DATA.usestdtemplates =0;
            DATA = AllV.CalcTemplateScores(DATA,fits{j}.Template);
        end
        sfits{j} = SlimFit(DATA, fits, j);
        nc(j,1) = sum(sfits{j}.SU > 0); 
        nc(j,2) = length(sfits{j}.SU);
    end
    details.nc = nc;
    DATA.autofits = sfits;
    if DATA.usestdtemplates == 0
        DATA.usestdtemplates =1;
        DATA = AllV.CalcTemplateScores(DATA);
    end
    [nmax,b] =  max(nc(:,1));
    if nmax ==2
%if "slimming" takes all fits down to 2 clusters, then
%use an unslimmed fit to start from
        [nmax,b] =  max(sum(nc,2));
    end

    sfits{b} = ConsensusFit(DATA, sfits{b},sfits,'ncrit',1);
    fit = sfits;
    return;
end    
%Under devel
% also try finding clusters common to everyone, then attaching the 
% ones left out.
if slimfits
    if slimfit > 0 || ~iscell(fits)
        fit = SlimFit(DATA,fits, slimfit);
        return;
    end
    for j = 1:length(fits)
        fit{j} = SlimFit(DATA, fits,j);
        sfits{j} = AllV.SetFit(fit{j},'slim');
        nc(j,1) = sum(sfits{j}.SU > 0); 
        nc(j,2) = length(sfits{j}.SU);
    end
    details.nc = nc;
    [~,b] =  max(nc(:,1));
    a = ConsensusFit(DATA, sfits{b},sfits,'ncrit',1);
    a.combinedfrom = b;
    fit = a;
    return;
end



for j = 1:length(fits)
    if isfield(fits{j},'combinedetails') %already combined
        nc(j) = 0;
    elseif isfield(fits{j},'SU')
        nc(j) = length(fits{j}.SU);
    else
        nc(j) = 1;
    end
end


[a,b] = max(nc);
minclusters = max([min(nc) 2]);
allgood = find(nc ==a);
checkfits = [];
if length(allgood) > 1
    fitscores = CellToMat(fits,'fitscore');
    if ~isempty(fitscores)
        [a,b] = max(fitscores(allgood,1));
        b = allgood(b);
    end

    for j = 1:length(allgood)
        for k = 1:j-1
            [matches, mdetails] = clust.MatchLsts(fits{allgood(j)}.clst,fits{allgood(k)}.clst);
            no = find(sum(mdetails.score > 0.8,2) >2);
            gmatch(j,k) = min(matches);
            if ~isempty(no)
                fprintf('Fits %d and %d have %d clusters, but different\n',allgood(j),allgood(k),a);
                checkfits(end+1,:) = [j k];
            end
        end
    end
end
maxpt = b;
details.maxfit = b;
fprintf('Max #Clusters in fit %d\n',b);
if consensus
    fit = ConsensusFit(DATA, fits{b},fits);
    return;
end


fit = fits{b};
pcspace = [];
tmpspace = [];
for j = 1:length(fit.SU)
        if fit.space{j}(1)  == 1
            pcspace = cat(2,pcspace,fit.space{j}(2:end));
        elseif fit.space{j}(1)  == 3
            tmpspace = cat(2,tmpspace,fit.space{j}(2:end));
        end
end

fisspace = GetFitSpace(DATA, fit);

clst = fit.clst;
finished = 0;
subs = [];
while finished == 0
    cls = unique(clst);
    for j = 1:length(cls);
        x = CalcIsolation(fitspace(:,fit.space{j,j}(2:end)),clst,cls(j));
        D(j,j) = x(1);
        nD(j,j) = x(3);
        for k =1:j-1;
            twospace = unique(cat(2,fit.space{j,j}(2:end),fit.space{j,k}(2:end),fit.space{k,k}(2:end)));
            x = CalcIsolation(fitspace(:,twospace),clst,cls(j), cls(k));
            D(j,k) = x(1);
            D(k,j) = x(2);
            nD(j,k) = x(3);
            nD(k,j) = x(4);
            if mean(x(1:2)) < 2.5 && mean(x(3:4)) < 2 
                cscore(j,k) = mean(x);
            end
%any time collapsing two produces one cluster that is
%better isolated by the "near" measure, this must be
% good - unless its a double trigger. Add CalcEfficacy to check
            cclst = ones(size(clst));
            cclst(clst == cls(j)) = 2;
            cclst(clst == cls(k)) = 2;
            cx = CalcIsolation(fitspace,cclst,2);
            if cx(3) > nD(j,j) && cx(3) > nD(k,k)
                cscore(j,k) = 5;                
            end
        end
    end
    D = fit.D(cls,cls);
    nd = size(D,1);
    xid = sub2ind(size(D),2:nd,2:nd);
    D(xid) = NaN;
    spaces = fit.space(cls,cls);
    [a,b] = min(D(2:end,2:end));
    [c,d] = min(nD(2:end,2:end));
    if size(D,1) <= minclusters || min(a) > 2
        finished = 1;
    else
        [~,c] = min(a);
        idx = [b(c)+1 c+1];
        fprintf('Combining %d into %d\n',cls(idx));
        clst(clst==cls(idx(1))) = cls(idx(2));
        subs(end+1,:) = idx(1:2);
        cscore(j,k) = 0;
        cscore(k,j) = 0;
        clear D;
        clear nD;
    end
end
details.subs = subs;
if ~isempty(subs) %were substitutiions
    newcls = unique(clst); %new list
    newid = find(ismember(cls,newcls));
    fit.isolation = fit.isolation(newid,:);
    fit.NDisolation = fit.NDisolation(newid,:);
    fit.bestscores = fit.bestscores(newid);
    fit.SU = fit.SU(newid);
    fit.space = spaces;
else
    fprintf('No clusters were combined\n');
end
id = find(D < 2.5);
nid = find(nD < 2.5);

separate = [];
for j = 1:length(fits)
    [matches, mdetails] = clust.MatchLsts(clst,fits{j}.clst);
    no = sum(mdetails.score > 0.8,2);
    id = find(no >1);
    for k = 1:length(id)
        matchid = find(mdetails.score(id(k),:) > 0.8);
        fprintf('Checking Fit %d cells %s against cell%d....',j,sprintf(' %d',matchid),cls(id(k)));
        d = [fits{j}.D(matchid(1),matchid(2)) fits{j}.D(matchid(2),matchid(1))];
        if min(d) > 2.5
            fprintf('Fit %d cl %d and %d separation %.2f,%.2f',j,matchid(1),matchid(2),d);
            separate(end+1,:) = [j matchid(1) matchid(2)];
        end
        fprintf('\n');
    end
    details.matchscore{j} = mdetails.score;
   if strcmp(plottype,'image')
       GetFigure('FitMatches');
       imagesc(mdetails.score);
       title(sprintf('Fit %d: %d clusters',j,nc(j)));
    end
end
cls = unique(clst);
id = find(diff(cls) > 1);
if ~isempty(id)
    j = id(1)+1;
    if sum(cls == j) == 0 %missing
        while j < max(cls)
            clst(clst == j+1) = j;
            j = j+1;
        end
    end
end
if strcmp(plottype,'image')
    j = DATA.cluster.bestfit(1);
    [matches, mdetails] = clust.MatchLsts(clst,fits{j}.clst);
    no = sum(mdetails.score > 0.8); %check if best fit tells us to combine split cells her
    id = find(no > 1);
    for j = 1:length(id)
        splits = find(mdetails.score(:,id(j)) > 0.8);
        for k = 1:length(splits)
            for c = 1:k-1
                d = [D(splits(k),splits(c)) D(splits(c),splits(k))];
                neard = [nD(splits(k),splits(c)) nD(splits(c),splits(k))];
                if min(d) < 2 || min(neard) < 1 %collapse these back
                    fprintf('Should coalesce %d and %d\n',splits(k),splits(c));
                end
            end
        end
    end
    GetFigure('FitMatches');
    imagesc(mdetails.score);
    title(sprintf('Fit %d: %d clusters',j,nc(j)));
end

details.separate = separate;
fit.clst = clst;
fit.D = D;
fit.nearD = nD;
fit.subs = details.subs;
fit.combinedetails = details;



function sfit =  SlimFit(DATA, fits, fiti, varargin)

verbose = 1;
if iscell(fits)
    fit = fits{fiti};
else
    fit = fits;
    if isfield(fit,'fitnumber')
        fiti = fit.fitnumber;
    end
end
sfit = fit;
if ~isfield(fit,'SU')
    sfit.SU = [];
end
if length(fit.D) < 3
    return;
end

pcspace = [];
tmpspace = [];
if fit.fitspace(1) == 1
    fitspace = DATA.pcs;
else
    fitspace = DATA.TemplateScores;
end

finished = 0;
substitutions = [];
clst = fit.clst;
newisolation = fit.isolation;
while finished == 0
    cmpisolation = [];
    joins = [];
    cls = unique(clst);
    a = [];
    neartest = [];
    cmpscore = ones([length(cls) length(cls) 2]) .*NaN;
    for j = 1:length(cls)
        for k = 1:j-1;
            sp = unique(cat(2,fit.space{j,j}(2:end),fit.space{k,k}(2:end)));
            fitspaces{j,k} = sp;
%Calculate isolation for the combined cluster, and for the
%individuals in the same space
            isolation = CalcIsolation(fitspace(:,sp),clst,cls([j k]));
            isolationj = CalcIsolation(fitspace(:,sp),clst,cls(j));
            isolationk = CalcIsolation(fitspace(:,sp),clst,cls(k));
            isolations(j,k,:) = isolation;
            cmpi = cat(2,isolationj([1 3]),isolationk([1 3]),isolation([1 3]));
            cmpisolation(j,k,:) = cmpi;
            a(j,k,1) =  cmpi(:,5) > cmpi(:,1) ;
            a(j,k,2) = cmpi(:,5) > cmpi(:,3);
            a(j,k,3) =  cmpi(:,6) > cmpi(:,2);
            a(j,k,4) = cmpi(:,6) > cmpi(:,4);
            neartest(j,k) = cmpi(6) - mean(cmpi([2 4]));
            cmpscore(j,k,1) = cmpi(5) - mean(cmpi([1 3]));
            cmpscore(j,k,2) = cmpi(6) - mean(cmpi([2 4]));
        end
    end
    [bi,ai] = find(squeeze(sum(a,3)) > 2 & neartest >0,1);
    if isempty(bi) || j < 3
        if finished == 0
            subs{1} = substitutions;
            finished = 1;
        end
        x = sum(cmpscore(2:end,2:end,:),3);
        [bi,ai] = find(x > -0.2);
        bi = bi+1;
        ai = ai+1;
    else
%if find more than one, just fixed the biggest gain first.        
        if length(bi) > 1
            x = sum(cmpscore,3);
            id = sub2ind(size(x),bi,ai);
            [~,maxi] = max(x(id));
            bi = bi(maxi);
            ai = ai(maxi);
        end
        I = [];
        if ai == 1 && size(a,3) > 2 %check if these are closer to a non-zero cluster
            for j = 1:size(a,1)
                I(j,:) = CalcIsolation(fitspace(:,sp),clst,cls(bi),cls(j));
            end
            I(bi,:) = Inf;
            [imin, id] = min(I(:,3));
            [imin, nid] = min(I(:,1));
            if id > 1
                if verbose
                    fprintf('Fit %d Setting %d->%d Rahter than 1\n',fiti,cls(bi),cls(id))
                end                
                ai = id;
                if id > bi
                    ai = bi;
                    bi = id;
                else
                    ai = id;
                end
            elseif nid > 1
                if verbose
                    fprintf('Fit %d Setting %d->%d Rahter than 1\n',fiti,cls(bi),cls(nid))
                end       
                if nid > bi
                    ai = bi;
                    bi = nid;
                else
                    ai = nid;
                end
            end
        end
        cmpi = squeeze(cmpisolation(bi,ai,:));
        isolation = squeeze(isolations(bi,ai,:));
        sp = fitspaces{bi,ai};
        d = squeeze(a(bi,ai,:));
        
        substitutions(end+1,:) = cls([bi ai]);
        cleared(cls(bi)) = cls(ai);
        clearscore(cls(bi),:) = squeeze(cmpisolation(bi,ai,:));
        if verbose
            fprintf('Fit %d Setting %d->%d D%.2f,%.2f->%.2f  nD %.2f,%.2f->%.2f\n',fiti,cls(bi),cls(ai),cmpi([1 2 5 3 4 6]));
        end
        clsa = unique(clst);
        clst(clst ==cls(bi)) = cls(ai);
        clsb = unique(clst);
        sid = find(ismember(clsa,clsb));
        newisolation = newisolation(sid,:);
        sid = find(clsb == cls(ai));
        newisolation(sid,:) = isolation(1:size(newisolation,2));
        clsid = find(ismember(cls,unique(clst)));
        cls = unique(clst);
        sfit.cleared = cleared;
    end
end
if ~isempty(substitutions)
    nalt = GetAltFit(sfit, 'slim');
    newspace = sfit.space;
    while(max(diff(cls)) > 1)
        for j = 2:length(cls)
            if cls(j)-cls(j-1) > 1 %gap
                substitutions(end+1,:) = [cls(j) -1];
                clst(clst==cls(j)) = cls(j)-1;
                newspace(cls(j)-1,:) = newspace(cls(j),:);
                newspace(:,cls(j)-1) = sfit.space(:,cls(j));
            end
        end
        cls = unique(clst);
    end
    sfit.alternates{nalt}.substitutions = substitutions;
    K = sfit.alternates{nalt};
    K.type = 'slim';
    K.isolation = newisolation;
    K.score = AllV.CalcFitScore(AllV.SetFit(fit,'slim'));
    sfit.alternates{nalt} = AddCalcs(DATA, K, clst, sfit,newspace);
end
nalt = GetAltFit(sfit, 'slim2');

if ~isempty(bi)
    replaced = [];
    sfit.alternates{nalt}.substitutions = cat(1,substitutions,cls([bi ai]));
    for j = 1:length(bi)
        cmpi = squeeze(cmpisolation(bi(j),ai(j),:));
        isolation = squeeze(isolations(bi(j),ai(j),:));
        clsa = unique(clst);
        if ~ismember(cls(ai(j)),clsa) %don't use a cell that is gone already
            if length(replaced) >= cls(ai(j)) && ismember(replaced(cls(ai(j))),clsa)
                ai(j) = find(cls == replaced(cls(ai(j))));
            end
        end
        if ismember(cls(ai(j)),clsa) %don't use a cell that is gone already
            clst(clst ==cls(bi(j))) = cls(ai(j));
            replaced(cls(bi(j))) = cls(ai(j));
            clsb = unique(clst);
            sid = find(ismember(clsa,clsb));
            newisolation = newisolation(sid,:);
            sid = find(clsb == cls(ai(j)));
            newisolation(sid,:) = isolation;
        else %should not get here
            x = find(sfit.alternates{nalt}.substitutions(:,1) == cls(bi(j)) & sfit.alternates{2}.substitutions(:,2) == cls(ai(j)));
            sfit.alternates{nalt}.err = sprintf('Cant replace %d with %d',cls(ai(j)),cls(bi(j)));
            cprintf('red', '%s\n',sfit.alternates{2}.err);
        end
    end
    sfit.alternates{nalt}.type = 'slim2';
    sfit.alternates{nalt}.isolation = newisolation;
    sfit.alternates{nalt}.score = AllV.CalcFitScore(AllV.SetFit(fit,'slim2'));
    sfit.alternates{nalt} = AddCalcs(DATA, sfit.alternates{nalt}, clst, sfit);
end


function newfit = AddCalcs(DATA, K, clst, fit, newspace)
    clsa = unique(clst);
    if nargin > 4
        fit.space = newspace;
    end
    fitspace = GetFitSpace(DATA, fit);
    for j = 1:length(clsa)
        xy = fitspace(:,fit.space{j,j}(2:3));
        [ell.xy, scores, details] = FindEllipse(xy,clst,'cluster',clsa(j)-1);
        ell = CopyFields(ell,details,{'scores' 'guess'});
        ell.allxy = mean(xy);
        K.ellipse(clsa(j)) = ell;
        dropi = AllV.CalcTriggerDrop(DATA, clst == clsa(j));
        K.dropi(clsa(j)) = dropi(3);
    end
    newfit = K;

function newfit = ConsensusFit(DATA, fit, fits, varargin)

ncrit = 1;

j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'ncrit')
        j = j+1;
        ncrit = varargin{j};
    end
    j = j+1;
end

for j = 1:length(fits)
    [matches, details] = clust.MatchLsts(fit.clst,fits{j}.clst);
    allmatches(j,1:length(matches)) = sum(details.score > 0.8);
    [a,b] = find(details.score > 0.8);
    clst = zeros(1,length(fit.clst));
    for k = 1:length(a)
        clst(details.matchid{a(k),b(k)}) = a(k);
    end
    clsts(j,:) = clst;
end
newfit = fit;
nc = max(fit.clst)+1;
nalt = GetAltFit(newfit, 'consensus');
clst = fit.clst;
missing = sum(clsts==0); %for each spike # of fits where it was not a match
if min(missing) >= ncrit
    ncrit = min(missing);
    fprintf('%sSetting Critierion for Consensus to %d missing fits\n',AllV.IDStr(DATA),ncrit);
end
clst(missing > ncrit) = nc;
if sum(missing > ncrit)
    newfit.space(nc,:) = fit.space(nc-1,:);
    newfit.space(1:size(fit.space,1),nc) = fit.space(:,nc-1);
    newfit.space(nc,nc) = fit.space(nc-1,nc-1);
end
if sum(clst ==1) == 0
    clst(clst == min(clst)) = 1;
end
%dont let loop below remove old members from isolation if a cluster
%gets moved completely
isolation = fit.isolation;
if sum(clst ==nc) %did change somethi
    newfit.alternates{nalt}.clst = clst;
    newfit.alternates{nalt}.type =  'consensus';
    cls = unique(clst);
    fitspace = GetFitSpace(DATA, fit,'exact');
    for j = 1:length(cls)
        isolation(j,:) = CalcIsolation(fitspace,clst,cls(j));
    end
    newfit.alternates{nalt}.isolation = isolation;
    newfit.alternates{nalt}.nocluster = nc;
    newfit.alternates{nalt} = AddCalcs(DATA, newfit.alternates{nalt}, clst, fit, newfit.space);
end

function nalt = GetAltFit(fit, type, varargin)

if isfield(fit,'alternates')
    nalt = 0;
    for j = 1:length(fit.alternates)
        if strcmp(fit.alternates{j}.type,type)
            nalt = j;
        end
    end
    if nalt == 0
        nalt = j+1;
    end
else
    nalt = 1;
end

function groups = GroupFits(fits, varargin)
method = 2;
for j = 1:length(fits)
    for k = 1:j-1;
        [matches, mdetails] = clust.MatchLsts(fits{j}.clst,fits{k}.clst);
        no = find(sum(mdetails.score > 0.8,2) >2);
        if diff(size(mdetails.score)) == 0
            gmatch(j,k) = min(matches);
        else
            gmatch(j,k) = 0;
        end           
        gmatch(k,j) = gmatch(j,k);
        gmatch(j,j) = NaN;
    end
end
if method == 1
    squishDistanceMatrix(gmatch);
else
    ng = 0;
    Gmatch = gmatch;
    while nansum(gmatch(:))
        [maxc,b] = max(gmatch(:));
        [a,b] = ind2sub(size(gmatch),b);
        gid = find(gmatch(a,:) >= 0.8);
        bgid = find(gmatch(:,b) > 0.8);
        if ~isempty(gid)
            ng = ng+1;
            gid = unique([bgid' gid]);
            groups(ng).group = gid;
            gscores = gmatch(gid,gid);
            groups(ng).scores = unique(sort(gscores(~isnan(gscores))));
            gmatch(gid,:) = 0;
            gmatch(:,gid) = 0;
        else
            gmatch(~isnan(gmatch)) = 0;
        end
    end
    [a,b] = find(isnan(gmatch)); %left overs with no grouping
    for j = 1:length(a)
        groups(end+1).group = a(j);
        groups(end).scores = 1;
    end
end
if method ==3 %both
    squishDistanceMatrix(gmatch);
end

function fitspace = GetFitSpace(DATA, fit, varargin)
%return space where fit was done, 
%GetFitSpace(DATA, fit, 'exact') gives only dims used
%Can get in trouble if some fits used stdTemplates 
%and others were Cluster derived - end up with wrong set of TemplateScores
exact = 0;

j = 1; 
while j <= length(varargin)
    if strncmpi(varargin{j},'exact',5)
        exact = 1;
    end
    j = j+1;
end

if fit.fitspace(1) == 1
    fitspace = DATA.pcs;
else
    fitspace = DATA.TemplateScores;
    if isfield(fit,'Template') && DATA.usestdtemplates 
%This fit used cell template, but current scores did not    
         T = AllV.mygetappdata(DATA,'AllTemplates');
         if isfield(T,'Scores')
             fitspace = T.Scores;
         end
    end
end
if exact
    dims = fit.fitspace(2:end);
    if max(dims) > size(fitspace,2)
        DATA = AddError(DATA,'-show','%s Fit %d Fitspace Dim %d not valid',AllV.IDStr(DATA),fit.fitnumber,max(dims));
        dims = dims(dims <= size(fitspace,2));
    end
    fitspace = fitspace(:,dims);
end
