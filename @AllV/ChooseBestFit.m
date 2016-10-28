function [bestfit, details] = ChooseBestFit(fits, varargin)
%[bestfit, details] = AllV.ChooseBestFit(fits, ....), Choose best of automatic fits;


type = 'maxisolation';
%type = maxcells uses min separation between cl and other cl (individually) as test
%'maxisolation' uses isolation(1) = cell vs rest

cellcrit = 2.5;
verbose =1;
DATA = [];
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'cluster') || isfield(varargin{j},'exptno')
        DATA = varargin{j};        
    elseif strcmp(varargin{j},'silent')
        verbose = 0;
    end
    j = j+1;
end

details.scores = [];
AllD = {};
groupid = [];
for j = 1:length(fits)
%only use isolation for cls = 2:end. isolation for cl1 means nothing, and
%can cause double counting....
    if isfield(fits{j},'bestscores')
        if strcmp(type,'maxcells')
            best(j,1:length(fits{j}.bestscores)) = sort(fits{j}.bestscores,'descend');
        else
            if isfield(fits{j},'isolation') && ~isempty(fits{j}.isolation)
                best(j,1:size(fits{j}.isolation,1)-1) = sort(fits{j}.isolation(2:end,1),'descend');
            else
                best(j,:) = 0;
            end
        end
%do same calculation for nearD measure.  size(nbest) needs to match size(best)
        if size(fits{j}.isolation,2) > 2
            nbest(j,1:size(fits{j}.isolation,1)-1) = sort(fits{j}.isolation(2:end,3),'descend');
        else
            nbest(j,1:size(best,2)) = 0;
        end
        AllD{j} = fits{j}.D;
        if isfield(fits{j},'SU')
            nclusters(j) = length(fits{j}.SU);
        else
            nclusters(j) = 0;
        end
    elseif isnumeric(fits{j}) %just givrn AllD
        AllD{j} = fits{j};
        nclusters(j) = size(fits{j},2);
    end
    if isfield(fits{j},'group') && ~isempty(fits{j}.group)
        groupid(j) = fits{j}.group(1);
    end
end

if isempty(best) %may only have one fit
    besti = find(nclusters > 1);
    if ~isempty(besti)
        bestfit = fits{besti(1)};
        bestfit.besti = besti(1);
        bestfit.bestfit = besti;
        bestfit.confidence =1;
    else
        bestfit = fits{1};
        bestfit.besti = 1;
        bestfit.bestfit = [1 NaN];
        bestfit.confidence = 0;
        details.err = sprintf('No Fits found clusters');
    end
    return;
end
%each row of best is isolation quality rank ordered for each fit
b = [];
score = zeros(1,size(best,1));
ascores{1} = zeros(1,size(best,1));
%make sure that at least one score is non-zero
lowcrit = 2;
if max(best(:)) < lowcrit
    lowcrit = max(best(:));
end
%for each cluster where one row has the best separation
%add len-column number. So doing best on #1 wins, but doint best
% on 2 and 3 is alwo goodl
for j = 1:size(best,2) %max # clusters
    [a(j), b(j)] = max(best(:,j)); %best for this cluster
    value = size(best,2); %column value
    if(a(j) > 2.5)
        bestid = find(best(:,j) > 2.5 & best(:,j) > a(j) * 0.99);
        score(bestid) = score(bestid) + value;
%also  add a reduced value to close seconds      
        id = find(best(:,j) > 2.5 & best(:,j) > a(j) * 0.9);
        id = setdiff(id,bestid);
        score(id) = score(id)+value/1.5;
    elseif(a(j) >= lowcrit) %give a little weight to good extra
        id = find(best(:,j) > a(j) * 0.99);
        score(id) = score(id) + value/10;
    end
    
    [a(j), c(j)] = max(nbest(:,j)); %best for this cluster
    value = size(nbest,2) - j;
    if(a(j) > 2.5)
        bestid = find(nbest(:,j) > a(j) * 0.99);
        ascores{1}(bestid) = ascores{1}(bestid) + value;
        id = find(nbest(:,j) > 2.5 & nbest(:,j) > a(j) * 0.9);
        id = setdiff(id,bestid);
        ascores{1}(id) = ascores{1}(id)+value/1.5;
    end
end

details.scores = score;
for j = 1:length(fits)
    [ascores{2}(j), x] = AllV.CalcFitScore(fits{j});
    for k = 2:length(x.scores)
        ascores{k+1}(j) = x.scores(k);
    end
end
details.ascores = ascores;
details.bestscores{1} = best;
details.bestscores{2} = nbest;

if isempty(b) %no fits worked
    bestfit = fits{1};
    bestfit.bestfit(1) = 1;
    bestfitl.besti = 1;
elseif length(b) ==1 || b(1) == b(2)
    bestfit = fits{b(1)};
    bestfit.bestfit(1) = b(1);
    bestfit.besti = b(1);
else
    [a,c] = max(score);
    bestfit = fits{c};
    bestfit.bestfit(1) = c;
    [aa,ac] = max(ascores{1});
    bestfit.bestfit(2) = ac;
    bestfit.besti = c;
end



for j = 2:length(ascores)
    fitscores(:,j-1) = ascores{j};
end
fitscores(:,end+1) = score;
fitscores(:,end+1) = ascores{1};
details.fitscores = fitscores;
%now rank order all the fits and see which is top most often
useid = intersect([1 3 4 5],find(sum(fitscores) > 0));
rankscores = fitscores(:,useid);
[~, ranki] = sort(rankscores);
ranks = repmat([1:size(rankscores,1)]',1,size(rankscores,2));
for j = 1:size(rankscores,1)
    npz(j) = sum(ranks(ranki==j));
end
zbest = max(npz);

details.npz = npz;
gid = find(npz == zbest); %can be ties
[a,c] = Counts(ranki(end,:));
gida = c(find(a == max(a)));
%find most frequent in either total rank (npz) or rank for each column
goodfits = cat(2,gid,gida);
[a,c] = Counts(goodfits);
gid = c(find(a == max(a))); %best fit/fits based on ranking
%gid is now most frequent cut and ties in rankings

details.zscore = sum(zscore(rankscores),2);
[z, gida] = max(details.zscore);
allgoodfits = ranki(end,:);
if z > 0    
    goodfits = cat(2, goodfits, gida);
    allgoodfits = cat(2, allgoodfits, gida);
else
    gida = gid(1);
end

if length(gid) == 1 && gid == gida  %simple case
    besti = gid;
else
    if sum(gid == gida)
        besti = [gida setdiff(gid,gida)];
    else
        besti = [gida gid];
        fprintf('No Clear Winner for fits out of %s',sprintf(' %d',besti));
        [m, md]= clust.MatchLsts(fits{besti(1)}.clst,fits{besti(2)}.clst);
        fprintf(': lowest Match %.2f\n',min(m));
        if diff(size(md.matchid)) == 0  && min(m) > 0.9
            fprintf('But the first two are equivalent\n');
        elseif length(besti) == 2 %mark uncertainty
            besti(3) = 0;
        end
    end    
end


%check that other fits really are differet
gi =  unique(goodfits);
gi = gi(gi ~= besti(1));
for j = 1:length(gi)
    [m, md]= clust.MatchLsts(fits{besti(1)}.clst,fits{gi(j)}.clst);
    fprintf('Fit %d lowest Match %.2f\n',j,min(m));
end
fitorder = [besti(besti>0) setdiff(unique(allgoodfits),besti)];
%record any fits that won on any scale, but not in besti
other = setdiff(1:length(details.zscore),fitorder);
[~,a] = sort(details.zscore(other),'descend');
details.fitorder = [fitorder other(a)];
bestfit = fits{besti(1)};
bestfit.besti = besti(1);
%bestfit.bestfit will be recorded in eckercluster, so can change fits
%without losing it..
if isfield(bestfit,'fitnumber')
    bestfit.bestfit(1) = bestfit.fitnumber;
    if isfield(bestfit,'alternateid')
        bestfit.bestfit(2) = bestfit.alternateid;
    else
        bestfit.bestfit(2) = 0;
    end
else
    bestfit.bestfit = besti(1);
end
bestfit.autofiti = bestfit.bestfit; %also the one uesd....
bestfit.eckercluster.AllD = AllD;
[bestfit.rankscore, details.rankorder] = Counts(allgoodfits);
if max(best(:)) == 0
    details.err = sprintf('No Fits  any good\n');
    bestfit.confidence = 0;
else
    bestfit.confidence = sum(allgoodfits == besti(1))./length(allgoodfits);
end
details.goodfits = setdiff(unique(allgoodfits),besti);

if ~isempty(groupid)
    details.alternatefit = [];
    g = groupid(details.fitorder(1));
    for j = 2:length(details.fitorder)
        fiti = details.fitorder(j);
        if ~ismember(groupid(fiti),g) && ascores{2}(fiti) > 4
            g = [g groupid(fiti)];
            if isfield(fits{fiti},'alternateid')
                details.alternatefit(end+1,2) = fits{fiti}.alternateid;
            else
                details.alternatefit(end+1,2) = 0;
            end
            details.alternatefit(end,1) = fits{fiti}.fitnumber;
            details.alternatefit(end,3) = fiti;
            details.alternatefit(end,4) = ascores{2}(fiti);
            bestfit.alternatefit = details.alternatefit;
        end
    end
end

a = unique(besti);
a = a(a>0);

if length(a) > 1
   [ details.matches, matchdetails] = clust.MatchLsts(fits{a(1)}.clst,fits{a(2)}.clst);
end

if verbose
    xstr = '';
    if isfield(bestfit,'alternateid')
        xstr = ['  ' bestfit.alternates{bestfit.alternateid}.type];
    end
    if isfield(bestfit,'fitnumber')
        fiti = bestfit.fitnumber;
    else
        fiti = bestfit.bestfit(1);
    end
    fprintf('%sBest Fit (%d%s) has %d clusters %d SU. Isolations:', AllV.IDStr(DATA),fiti,xstr,length(bestfit.SU),sum(bestfit.SU));
for j = 1:length(bestfit.SU)
    if bestfit.SU(j)
        fprintf(' %.3f (%.3f)',bestfit.isolation(j,1),bestfit.nearD(j,1));
    else
        fprintf(' *%.3f (%.3f)',bestfit.isolation(j,1),bestfit.nearD(j,1));
    end
end
fprintf('\n');
end
