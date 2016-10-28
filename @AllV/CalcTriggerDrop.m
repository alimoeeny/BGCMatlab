function [dropi, gfit, details] = CalcTriggerDrop(DATA, clid, varargin)
%CalcTriggerDrop(DATA, clid) estimate Dropped spikes from trigger hist
%Also works if DATA is a Cluster struct containing triggerV and clst
if length(clid) > 1
    id = clid;
    nid = setdiff(1:length(DATA.clst),id);
else
    id = find(DATA.clst == clid+1);
    nid = find(DATA.clst ~= clid+1);
end
ncut = length(id);
details.vhist = 0;

if isfield(DATA,'triggerV') % A cluster struct
    rV = DATA.triggerV;
    DATA.interactive = 0;
else
    rV = DATA.rV;
end
if cellstrcmp('noplot',varargin)
    DATA.interactive = -1;
end
if length(id) < 10
    nbins(1) = 1;
    d = (mean(rV(id)) - DATA.Trigger(1))./std(rV(id));
    if DATA.Trigger(1) < 0
        d = -d;
    end
    if d^3 < 2.2 %adjust for bias at low values. see sims/dropi
        d = 0;
    else
        d = (d.^3 - 2.3).^(1/3);
    end
    dropi = [0 0 d NaN];
elseif length(id) < 100
    nbins(1) = 10;
elseif length(id) < 5000
    nbins(1) = round(length(id)./20);
else
    nbins(1) = 250;
end
if length(nid) < 10
    nbins(2) = 1;
elseif length(nid) < 100
    nbins(2) = 10;
elseif length(nid) < 5000
    nbins(2) = round(length(nid)./20);
else
    nbins(2) = 250;
end



if nbins(1) > 2 && ~isempty(id)
    V = rV(id);
    [a,b] = hist(V,nbins(1));
    cluster.vhist = a;
    cluster.vhistrange = minmax(b);
    [c,d] = hist(rV(nid),nbins(2));
    details.vhist = [a(:) b(:)];
    details.muvhist = [c(:) d(:)];
    cmax = max(a);
    
    if DATA.Trigger(1) < 0
        tid = find(b <= DATA.Trigger(1));
        ntid = find(d <= DATA.Trigger(1));
    else
        tid = find(b >=DATA.Trigger(1));
        ntid = find(d >= DATA.Trigger(1));
    end
    
    if length(tid) < 2
        dropi = [0 0 0 NaN];
        gfit.sd = NaN;
    else
        gfit = FitGauss(b(tid),a(tid),'maxiter',500);
        if nbins(2) > 1
            ngfit = FitGauss(d(ntid),c(ntid),'maxiter',500);
        else
            ngfit.sd = 0;
            ngfit.mean = d;
        end
        if length(tid) && DATA.interactive >= 0
            hold on;
            plot(b(tid),gfit.fitted,'k');
        end
        nclose = sum(abs(V-DATA.Trigger(1)) < std(V)/10);
        dropi(1) = nclose./ncut;
        dropi(2) = abs(mean(V)-DATA.Trigger(1))./std(V);
        if DATA.Trigger(1) < 0
            dropi(3) = (DATA.Trigger(1)-gfit.mean)./abs(gfit.sd);
            dropi(4) = (DATA.Trigger(1)-ngfit.mean)./abs(ngfit.sd);
        else
            dropi(3) = (gfit.mean-DATA.Trigger(1))./abs(gfit.sd);
            dropi(4) = (DATA.Trigger(1)-ngfit.mean)./abs(ngfit.sd);
        end
    end
else
    dropi = [0 0 0 NaN];
    gfit.sd = NaN;
end