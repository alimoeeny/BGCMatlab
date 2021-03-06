function [cp, details] = GrandCP(Expt, nmin, varargin)

exa = Expt.Stimvals.et;
exb = Expt.Stimvals.e2;
twin = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'times',5)
        j = j+1;
        twin = varargin{j};
    end
    j = j+1;
end

id = find(abs([Expt.Trials.RespDir]) ==1);
Trials = Expt.Trials(id);
avals = [Trials.(exa)];
bvals = [Trials.(exb)];
choices =  [Trials.RespDir];
av = unique(avals);
bv = unique(bvals);
if length(twin) ==2
    for j = 1:length(Trials)
        Trials(j).count = sum(Trials(j).Spikes >= twin(1) & Trials(j).Spikes <= twin(2));
    end
end
counts = [Trials.count];

nx = 0;
for j = 1:length(av)
    for k = 1:length(bv)
        if abs(av(j)) > 0.0001
            if k == 1
            id = find(avals ==av(j));
            else
                id = [];
            end
        else
            id = find(avals ==av(j) & bvals == bv(k));
        end
        if sum(choices(id) == -1) > nmin && sum(choices(id) == 1) >= nmin;
            nx = nx+1;
            zscores{nx} = (counts(id)-mean(counts(id)))./std(counts(id));
            zchoices{nx} = (choices(id) ==1); %0 or 1
            details.zcp(nx) = CalcCP(counts(id(choices(id) ==1)),counts(id(choices(id) ==-1)));
        end
        ns(j,k) = length(id);
    end
end

pcounts = [];
ncounts = [];
for j = 1:nx
    pcounts = [pcounts, zscores{j}(zchoices{j})];
    ncounts = [ncounts, zscores{j}(~zchoices{j})];
end
cp = CalcCP(pcounts,ncounts);

nresample = 1000;
for n = 1:nresample
    pcounts = [];
    ncounts = [];
    for j = 1:nx
        C = zchoices{j}(randperm(length(zchoices{j})));
        pcounts = [pcounts, zscores{j}(C)];
        ncounts = [ncounts, zscores{j}(~C)];
    end
    cps(n) = CalcCP(pcounts,ncounts);
end
details.pval = 2 * min([length(find(cps > cp)) length(find(cps < cp))])/(nresample);
details.cps = cps;
details.ns = ns;
details.choicesign = mean(choices(avals>0))-mean(choices(avals < 0))
