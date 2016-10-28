function [DATA, C, basefit] = SetClusterFromFit(DATA, C, fit, varargin)
%[DATA, C, basefit] = AllV.SetClusterFromFit(DATA, C, fit, 
%copies isolation, dropi, from fit to Cluster C
refitellipse = 0;
state.verbose = 0;
showquality = 0;
clusterprops = [];
if DATA.interactive > 0
    state.verbose = 1;
end

if isnumeric(fit)
    fit = AllV.SetFit(C.autofits, fit);
end
basefit = fit; 



j =1;
while j <= length(varargin)
    if isfield(varargin{j},'combinecls')
        clusterprops = varargin{j};
    elseif strncmpi(varargin{j},'nofit',5)
        refitellipse = 0;
    elseif strncmpi(varargin{j},'quiet',5)
        state.verbose = 0;
    elseif strncmpi(varargin{j},'refit',5)
        refitellipse = 1;
    elseif strncmpi(varargin{j},'showquality',5)
        showquality = 1;
    elseif sum(strcmp(varargin{j},{'slim' 'slim2' 'consensus'}))
        fit = AllV.SetFit(fit, varargin{j});
    end
    j = j+1;
end
if state.verbose
    fprintf('Setting parameters in Cluster from Kalman Fits...');
end
C = DATA.cluster;
if isfield(fit,'eckercluster')
    C.eckercluster = fit.eckercluster;
end
C.eckercluster = CopyFields(C.eckercluster,fit,{'energy' 'fitspace' 'bestfit' 'alternatefit' 'SU'});
C = CopyFields(C, fit, {'bestfit' 'autofiti'});
C.next = {};  %clear this out. Will be rebuilt below
C.ctime = now;
if isfield(C,'bestfit') && length(C.bestfit) == 1
    C.bestfit(2) = 0;
end
if isfield(clusterprops,'combinecls')
    cmb = clusterprops.combinecls+1;
    fit.clst(fit.clst == cmb(2)) = cmb(1);
    if length(C.eckercluster.SU) >= cmb(2)
        C.eckercluster.SU(cmb(2)) = [];
    end
    if length(C.eckercluster.sortorder) >= cmb(2)        
        C.eckercluster.sortorder(cmb(2)) = [];
        C.eckercluster.type1(cmb(2)) = [];
        C.eckercluster.type2(cmb(2)) = [];
    else
        fprintf('Requested combine of non-existent clusters %s %d,%d\n',str(C),cmb);
    end
    C.eckercluster.combinecls = cmb; 
    if isfield(C,'bestfit')
        C.eckercluster.autochoice = C.bestfit;
        C.autochoicefit = C.bestfit;
    end
end

DATA.clst = fit.clst;
if ~isfield(fit,'SU') %failed
    fit.SU = 0;
end
mu = find(fit.SU ==0);

DATA.clid = find(DATA.clst == C.cluster+1);
DATA.nid = find(DATA.clst ~= C.cluster+1);
DATA.gmcid = DATA.clst;
cls = unique(DATA.clst);
if isfield(DATA,'t')
    clt = DATA.t;
elseif isfield(DATA.cluster,'t')
    clt = DATA.t;
end
spaces = [];

Tbins = linspace(min(DATA.rV),max(DATA.rV));

if DATA.interactive > 0
    fargs = {'plot' 'tag' 'FindEllipse'};
else
    fargs = {};
end
%clst = 1 should always be MU. so only do this for > 1
cls = cls(cls>1);
DATA.xy = {};
%a == a is MU, so only cls-1 clusters
if isfield(fit,'Template') && DATA.usestdtemplates == 1
    DATA.cluster.neednewtemplate = 1;
    DATA.usestdtemplates = 0;
    DATA = AllV.CalcTemplateScores(DATA, fit.Template);
elseif ~isfield(fit,'Template') && DATA.usestdtemplates == 0
    DATA.usestdtemplates = 1;
    DATA = AllV.CalcTemplateScores(DATA);
end

%This is just in case. xy{1} shoudl be set byt fit{2,2} (clst ==2 is
%cluster 1
if ~isempty(cls)
    fit.bestspace(1,:) = fit.space{1,1}(1:3);
    if fit.bestspace(1,1) == 1
        DATA.xy{1} = DATA.pcs(:,fit.bestspace(1,2:3));
    elseif fit.bestspace(1,1) == 3
        DATA.xy{1} = DATA.TemplateScores(:,fit.bestspace(1,2:3));
    end
end



%Mapping from fit.clst to cluster is truly arbitray
% could even consider re-ordering at this point.  But in any case skip
% missing entries. = use j not c, to set cluster element.
%i.e if there are no clst ==2, cluster.next{1} is clst ==4
othermeans = {};
for j = 1:length(cls)   
    %diagonal in fit.space is cluster j vs all others
    c = cls(j); %element of clst;
    spaces = [spaces fit.space{j+1,j+1}(1:3)];
    fit.bestspace(j+1,:) = fit.space{j+1,j+1}(1:3);
    if fit.bestspace(j+1,1) == 1
        DATA.xy{j} = DATA.pcs(:,fit.bestspace(j+1,2:3));
    elseif fit.bestspace(j+1,1) == 3
        DATA.xy{j} = DATA.TemplateScores(:,fit.bestspace(j+1,2:3));
    end
    sid = find(DATA.clst == c);
    if c > j+1 %skipped one
        DATA.clst(DATA.clst == c) = j+1;
    end
    clst = ones(size(DATA.clst));
    clst(sid) = 2;
    if isempty(sid) %no spikes in this group - should not happen
        cprintf('red','No spikes for clst = %d\n',c); 
        fit.ellipse(j).xy = [NaN NaN NaN NaN NaN];
    elseif ~isfield(fit,'ellipse') || length(fit.ellipse) < j+1 || refitellipse
        [ell.xy, scores, details] = FindEllipse(DATA.xy{j},clst,fargs{:});
        ell = CopyFields(ell,details,{'scores' 'guess'});
        ell.allxy = mean(DATA.xy{j});
        fit.ellipse(j+1) = ell;
    else
        ell = fit.ellipse(j+1);
    end
    [a,b] = hist(DATA.rV(sid),Tbins);
%if clst has gaps, close these up here. Just need to watch calculation of
%dropi and isolation later
    if j > 1 
        C.next{j-1}.space = fit.bestspace(j+1,:);
        C.next{j-1}.MeanSpike = AllV.PlotMeanSpike(DATA,'recalc','cluster',j,'noplot');
        if isfield(fit,'eckercluster')
            C.next{j-1}.mahal(1) = fit.eckercluster.type1(j+1);
            C.next{j-1}.mahal(2) = fit.eckercluster.type2(j+1);
            C.next{j-1}.mahal(3) = fit.eckercluster.SU(j+1);
            C.next{j-1}.isSU = fit.eckercluster.SU(j+1);
        else
            C.next{j-1}.mahal(3) = fit.SU(j+1);
            C.next{j-1}.isSU = fit.SU(j+1);
        end
        C.next{j-1}.mahal(4) = 0;
        C.next{j-1}.vhist = a;
        C.next{j-1}.auto = 3;
        C.next{j-1}.shape = 3;
        C.next{j-1}.space = fit.space{j+1,j+1};
        C.next{j-1}.isolation = fit.isolation(j+1,:);       
        C.next{j-1}.bestisolation.isolation = fit.isolation(j+1,:);       
        C.next{j-1}.bestisolation.space = fit.bestspace(j+1,:);
        C.next{j-1}.NDisolation = fit.NDisolation(j+1,:);
        C.next{j-1}.vhistrange = minmax(b);
        C.next{j-1}.cluster = j; %in case clst has gaps
        C.next{j-1}.ellipse = ell;
        if isfield(fit,'optellipse')
            C.next{j-1}.optellipse = fit.optellipse(j+1);
        end
        C.next{j-1}.ctime = now;
        C.next{j-1}.times = clt(sid);
        C.next{j-1}.angle = 0;
        C.next{j-1}.dropi = AllV.CalcTriggerDrop(DATA, cls(j)-1);
        fit.dropi(c) = C.next{j-1}.dropi(3);
        Amps(j+1,:) = C.next{j-1}.MeanSpike.amp./max(C.next{j-1}.MeanSpike.amp);
        othermeans{j-1} = C.next{j-1}.MeanSpike.ms;
    else
        C.MeanSpike = AllV.PlotMeanSpike(DATA,'recalc','cluster',1,'noplot');
%need this line otherwise cluster.space isn't correct. E.G. when loading pprM054 E53)3        
        C.space = fit.bestspace(2,:);
        C.vhist = a;
        C.vhistrange = minmax(b);
        C.isolation = fit.isolation(2,:);       
        C.bestisolation.isolation = fit.isolation(2,:);       
        C.bestisolation.space = fit.bestspace(2,:);
        if isfield(fit,'NDisolation')
            C.NDisolation = fit.NDisolation(2,:);
        end
        C.ellipse = ell;
        if isfield(fit,'optellipse')
            C.optellipse = fit.optellipse(2);
        end
        if isfield(fit,'eckercluster')
            C.isSU = fit.eckercluster.SU(j+1);
        else
            C.isSU = fit.SU(j+1);
        end
        Amps(2,:) = C.MeanSpike.amp./max(C.MeanSpike.amp);
        C.angle = 0;
        C.times = clt(sid);
        C.dropi = AllV.CalcTriggerDrop(DATA, 1);
        fit.dropi(c) = C.dropi(3);
        C.cluster = c-1;
    end        
end

%if we have re-calculated ellipses for an altenate fit
% record this
if isfield(fit,'ellipse')
    if ~isfield(fit,'dropi')
        fit.dropi = NaN;
    end
    if isfield(fit,'ellipse') && isempty(fit.ellipse(1).xy)
        fit.ellipse(1).xy = [NaN NaN NaN NaN];
    end
    if isfield(fit,'alternateid')
        basefit.alternates{fit.alternateid}.ellipse = fit.ellipse;
        basefit.alternates{fit.alternateid}.dropi = fit.dropi;
    else
        basefit.ellipse = fit.ellipse;
        basefit.dropi = fit.dropi;
    end
end
if showquality
    fprintf('\n');
    [ns, cls] = Counts(fit.clst);
    if ~isfield(fit,'energy')
        fit.energy = zeros(size(cls));
    end
    for j = 1:min([size(fit.isolation,1) length(fit.SU)]);
        if fit.SU(j)
            su = '*';
        else
            su = '';
        end
        fprintf('%d%s: %.2f %.2f dropi %.2f En%.1f, N%d\n',cls(j),su,fit.isolation(j,[1 3]),fit.dropi(j),fit.energy(j),ns(j));
    end
    [a,b] = AllV.CalcFitScore(fit);
    fprintf('Fit Score:')
    for j = 1:length(b.scores)
        fprintf(' %.2f',b.scores(j));
    end
    fprintf('\n');
    if isfield(fit,'errs')
        for j = 1:length(fit.errs)
            cprintf('red','%s\n',fit.errs{j});
        end
    end
end

if length(cls) == 1
    C.angle = 0;
    C.isSU = 0;
end
sid = find(DATA.clst == 1);
if ~isfield(fit,'space')
    fit.bestspace = [1 1 2];
    fit.isolation = 0;
elseif length(cls) <= 1 %no cluster found
    spaces = [spaces fit.space{1,1}(1:3)];
    fit.bestspace(1,:) = fit.space{1,1}(1:3);
    if fit.bestspace(1,1) == 1
        DATA.xy{1} = DATA.pcs(:,fit.bestspace(2:3));
    elseif fit.bestspace(1,1) == 3
        DATA.xy{1} = DATA.TemplateScores(:,fit.bestspace(1,2:3));
    end
    C.dropi = [0 0 NaN NaN];
    fit.bestspace(1,:) = fit.space{1,1}(1:3); %best space for 0 vs 1
else
    fit.bestspace(1,:) = fit.space{2,2}(1:3); %best space for 1 vs all
    for j = 2:size(fit.space,1)
        fit.bestspace(j,:) = fit.space{j,j}(1:3);
    end
end

V = AllV.mygetappdata(DATA,'AllVoltages');
ms = mean(V(:,:,sid),3);
A = std(ms');
Amps(1,:) = A/max(A);
C.eckercluster.bestspace = fit.bestspace;
C.eckercluster.isolation = fit.isolation;
C.eckercluster.amplitudes = Amps;

C.autocutmode = 'ecker';
C = CopyFields(C,DATA,{'probe','Trigger' 'autocutmode' 'chspk'});
C.auto  = 3;
C.MeanSpike = AllV.PlotMeanSpike(DATA,'recalc','noplot');
if ~isempty(othermeans)
    good = [];
    for j = 1:length(othermeans)
        if ~isempty(othermeans{j})
            good(j) = 1;
        else
            good(j) = 0;
        end
    end
    othermeans = othermeans(find(good));
    C.MeanSpike.othermeans = othermeans;
    for j = 1:length(othermeans)
        xc = corrcoef(C.MeanSpike.ms,othermeans{j});
        C.MeanSpike.otherxc(j,1) = xc(1,2);
        xc = corrcoef(C.MeanSpike.mu,othermeans{j});
        C.MeanSpike.otherxc(j,2) = xc(1,2);
    end
end
C.clst = DATA.clst;
C.bmc = 0;
if isfield(fit,'eckercluster')
    C.mahal(1) = fit.eckercluster.type1(1);
    C.mahal(2) = fit.eckercluster.type2(1);
    C.mahal(3) = fit.eckercluster.SU(1);
else
    if ~isempty(fit.SU)
        C.mahal(3) = fit.SU(1);
    else
        C.mahal(3) = 0;
    end        
end
C.mahal(4) = 0;
C.eDistance = fit.D;
C.ctime = now;
C.shape = 3;

DATA.usegmcid = 1;  %plot results according to GM clustering
C = CopyFields(C,DATA,{'TemplateUsed' 'TemplateB' 'mumeanUsed' 'DprimeUsed' 'usestdtemplates'},'-noempty');


details.newDATA = 0;
C.cluster = 1;
C.ncut = length(DATA.clid);
C.nspks = length(DATA.clst);
C = rmfields(C,'jamescluster');
C.autocutmode = 'ecker';
if ~isfield(C,'autofiti') || length(C.autofiti) == 1
    C.autofiti = [0 1];
end
C.ctime = now;

if DATA.interactive > 0
    fprintf('\n');
end
