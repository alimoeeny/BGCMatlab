function [shakes, details] = FindShake(FX, varargin)
%shakes = PC.FindShakes(AutoList) use rate jumps to find
%shakes, so that don't waste time looking at cluster alignment,

if isfigure(FX)
    F = FX;
    FX = getappdata(1,'AutoCellList');
    R = FX.ratecheck;
elseif isfield(FX, 'blkcv')
    R = FX;
else
    R = FX.ratecheck;
end
expts = [];
sizes = [];
jumpids = [];
%Expt type, and cell no
[e,p] = find(R.njump > 0);
for j = 1:length(e)
    jx = R.jumps(e(j),p(j));
    expts = cat(2, expts, jx.ex);
    sizes = cat(2, sizes, jx.sz);
    jumpids = cat(2, jumpids, repmat([e(j) p(j)]',1,size(jx.ex,2)));
end
ncrit = length(R.cellid)./4;
if ncrit < 3
    ncrit = 3;
end
[e,p] = find(abs(R.step) > 0);
stepids = [];
stepex = [];

for j = 1:length(e)
    jx = R.steps{e(j),p(j)};
    [stepsz(j), id] = max(abs([jx.step]));
    stept(j) = jx(id).t;
    stepex(j) = jx(id).e;
    stepids(:,j) = [e(j) p(j)];
end

[njs,b] = Counts(expts);
njexpts = b;
for j = 1:length(b)
    sz(j) = sum(abs(sizes(expts(2,:)== b(j) | expts(2,:) ==b(j))));
    nj(j) = sum(expts(1,:)== b(j) | expts(2,:) ==b(j));
end
plot(b,sz,'o');


id = find(stepsz > 0.4);
[nstep,b] = Counts(stepex(id));
stepexpts = b;
tq(1:max(b)) = NaN;
stepdata = [];
nx = 0;
for j = 1:length(b)
    xid = find(stepex(id) ==b(j));
    t = stept(id(xid));
    tsd(j) = std(t);
    cid = id(xid);
    for c = cid(:)'
        e = stepids(1,c);
        p = stepids(2,c);
        jx = R.steps{stepids(1,c),stepids(2,c)};
        sid = find(abs([jx.step]) > 0.4);
        for k = 1:length(sid)
            stepdata(:,end+1) = [e p sid(k) jx(sid(k)).step jx(sid(k)).e jx(sid(k)).t];
        end
    end
    nx = nx+length(xid);
    nxs(j) = length(xid);
end
for j = 1:length(b)
    xid = find(stepdata(5,:) == b(j));
    if length(xid) >= ncrit
        tq(b(j)) = diff(prctile(stepdata(6,xid),[25 75]));
    end
end
badj = find(tq < 5); %ncrit cells with jumps within 5 trials
badstep = b(find(nstep >= ncrit));
bade = njexpts(find(njs > ncrit));
badexpts = union(badj, badstep);
badexpts = union(badexpts, bade);
ns = 0;
for j = 1:length(badexpts)
    e = badexpts(j);
    aid = find(njexpts == e);
    bid = find(stepexpts == e);
    fprintf('E%d %d jumps, ',badexpts(j),njs(aid));
    xid = find(expts(1,:)== e | expts(2,:) == e);
    n = [0 0 0];
    nsm =  0;
    for k = 1:length(xid)
        pos = jumpids(:,xid(k));
        jx = R.jumps(pos(1),pos(2));
        n(jx.type+1) = n(jx.type+1) +1;
    end
    fprintf(' (%dff %dstep %d) ',n(2),n(3),n(1));
    fprintf('%d steps',nstep(bid));

    
    xid = find(stepdata(5,:) == e);
    if ~isempty(xid)
        probes = stepdata(2,xid);
        nprobes = length(unique(probes));
        shaket = ceil(prctile(stepdata(4,xid),50));
        fprintf(' sz%.2f %d probes tscatter %.1f (%d)',nanmean(abs(stepdata(4,xid))),nprobes,tq(e),shaket);
    else
        nprobes = 0;
    end
    fprintf('\n');
    if isfield(FX,'toplevel') && isfigure(FX.toplevel);
        F = GetFigure(sprintf('Expt%dRates',e),'parent',FX.toplevel);
        PC.PlotClusterRates(FX.toplevel,'rateseq','expt',e,'overlap',2,'figure',F,'quiet');
    end
end
    details(j).nexpt = njs(aid);
    details(j).stepdata = stepdata;
    details(j).expt = e;
    if njs(aid) > 10 && tq(e) < 10  && nprobes > ncrit
        ns = ns+1;
        shakes(ns).expt = e;
        shakes(ns).trial = shaket;
    end
end