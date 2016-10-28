function [xc, details] = xcorrtimes(at, bt, varargin)
%[xc, details] = xcorrtimes(a, b, varargin)
%return the cross-correlogram between two lists of timestamps (a, b)  
%in second. By default binwith is 1ms, range is +- 200ms
%xcorrtimes(a,b, 'times', t) sets vector of times for histogram

times = -0.2005:0.001:0.2005;
btimes = [];
clipping = 0;
method = 7; %was 2 before
maxt = max(abs(times));
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'clip',4)
        clipping = 1;
    elseif strncmpi(varargin{j},'btimes',4)
        j = j+1;
        btimes = varargin{j};
    elseif strncmpi(varargin{j},'maxt',4)
        j = j+1;
        maxt = varargin{j};
    elseif strncmpi(varargin{j},'method',4)
        j = j+1;
        method = varargin{j};
    elseif strncmpi(varargin{j},'selftest',4)
        ns = [1000 10000];
        for j = 1:length(ns)
            n = ns(j);
            x = round(rand(n,2) * 1000);
            tic;
            xc = xcorrtimes(x(:,1),x(:,2));
            ts(j,1) = toc;
            tic;
            xc = xcorrtimes(x(:,1),x(:,2),'method',2);
            ts(j,2) = toc;
            tic;
            xc = xcorrtimes(x(:,1),x(:,2),'method',3);
            ts(j,3) = toc;
        end
        details.ts = ts;
        details.ns = ns;
    elseif strncmpi(varargin{j},'times',4)
        j = j+1;
        times = varargin{j};
    end
    j = j+1;
end

details.counts = [length(at) length(bt)];
details.xpts = (times(1:end-1)+times(2:end))/2;
times = reshape(times,1,length(times));
xc = zeros(size(times));
if method == 2    
    bt = reshape(bt,1,length(bt));
    for j = 1:length(at);
        dt =  at(j) - bt;
        cc = histc(dt(abs(dt)<= maxt),times);
        if ~isempty(cc)
            xc = xc + cc;
        end
    end
    [a,b] = min(abs(times));
    details.efficacy = [xc(b)./length(at) xc(b)./length(bt)];
    details.midpt = b;
    xc = xc(1:end-1);
elseif method == 3
    for j = 1:length(at);
        dt(j,:) =  at(j) - bt;
    end
    xc = hist(dt(:),times);
elseif method == 4  % tests
    dts = [];
    ts = [];
    xc = [];
    for j = 1:length(at);
        dt =  at(j) - bt;
        id = find(abs(dt) < maxt);
        if length(id)
            dts = [dts dt(id)];
            ts = [ts bt(id)];
        end
    end
    xc(1,:) = dts;
    xc(2,:) = ts;
elseif method == 5
    id = find(diff(bt) > maxt);
    if isempty(bt)
        return;
    end
    if isempty(id)
        id(1) = length(bt);
    else
        id(end) = length(bt);
    end
    ct = bt(1:id(1));
    n = 1;
    at = at(at <= bt(end)+maxt);
    for t = at
        if t > ct(end)+maxt
            ct  = bt(id(n):id(n+1));
            n = n + 1;
        end
        dt =  t - ct;
        cc = histc(dt(abs(dt)< maxt),times);
        if length(cc)
        xc = xc + cc;
        end
    end
elseif method == 6 %may be faster - check with long vectors
    a = bsxfun(@minus,at,bt');
    xc = histc(a(:),times);
elseif method == 7  %split into semgents where gaps are large betweem
%    [a,b,c] = mymemory;
    ns = 1;
    maxsz = 1e6;
%first limit to overlapping portion of x,y    
    xc = zeros(size(times,1),size(times,2)-1);
    [a,b] = min(abs(times));
    details.midpt = b;
    details.efficacy = [0 0];
    if isempty(bt)
        return;
    end
    at = at(at> min(bt)-10 & at < max(bt)+10);
    if isempty(at)
        return;
    end
    bt = bt(bt> min(at)-10 & bt < max(at)+10);
    splitpt = clust.SplitTimes(at,bt);
    for j = 1:size(splitpt,1)-1
        xi = at(splitpt(j,1):splitpt(j+1,1));
        yi = bt(splitpt(j,2):splitpt(j+1,2));
        a = bsxfun(@minus,xi,yi');
%        xid = find(abs(a) < 0.001);
%histc(a(:)') gives the same size(1) even if length xi/yi == 1
        a = a(:)';
        as{j} = histc(a,times);
        if ~isempty(btimes)
            bs{j} = histc(a,btimes);
        end
    end
    xc = sum(cat(1,as{:}),1);
    xc = xc(1:end-1);
    b = details.midpt;
    adjust = (sum(xc) - xc(b))./(length(xc)-1);
    details.efficacy = (xc(b)-adjust)./[length(at)-adjust length(bt)-adjust];

    if ~isempty(btimes)
        bxc = sum(cat(1,bs{:}));
        bxc = bxc(1:end-1);
        details.btimes = (btimes(1:end-1)+btimes(2:end))/2;
        details.xcb = bxc;
    end
else
    for j = 1:length(at);
        dt =  at(j) - bt;
        cc = hist(dt,times);
        xc = xc + cc;
    end
end
if clipping
    xc = xc(2:end-1);
    details.xpts = details.xpts(2:end-1);
end
