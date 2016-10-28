function [e, details] = CalcEfficacy(x,y,varargin)
% [e, details] = CalcEfficacy(x,y,varargin)
%find peak in xcorrelogram over +- 1ms range.
%then count sync spikes +- 0.2ms from this
%e(1) is %of spikes in x within 0.2ms of an event in y = efficacy x-> y
%e(2) is %of spikes in y within 0.2ms of x
%e(1:2) con-incidence #s are adjusted by baseline rate.
%e(3) and e(4) are the same calculation for spikes <1ms apart but more than
%0.2ms. So if e(3) is as big as e(1) its not true synchrony
% negative e(3) or e(4) means spikes in the central ms outside 0.2 of peak
%are less common than expected from mean rate.
% e(5) and (6) are adjusted not from mean rate but from remaining counts in
%central ms.  So effects due to long term rate changes are discounted.
%so large values of e(5) and e(6) at least reflect the size of the peak in
%the central 0.2 ms
%x,y can be lists of itmes, or cluster structs (or any struct with field 'times')
details = [];
method = 2;

splits = 100;
ts = [-0.01005:0.0001:0.01005];

if isfield(x,'times')
    x = x.times;
end

if isfield(y,'times')
    y = y.times;
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'method',5)
        j = j+1;
        method = varargin{j};
    elseif strncmpi(varargin{j},'nosplit',5)
        splits = 1;
    end
    j = j+1;
end
e = [0 0 0 0 0 0 0 0];
details.delay = NaN;
if length(y) < splits * 2;
    splits = 1;
end
isi(1) = median(diff(x));
isi(2) = median(diff(y));
p = 1./(isi * 1000); %p(spike) in 1ms bin
nspk = [length(x) length(y)];
adjust(1)  = length(x) .* p(2);
adjust(2)  = length(y) .* p(1);
details.nspk = nspk;
if isempty(x) || isempty(y)
    return;
end


bindur = 2/10;
if method == 2
    ns = 1;
%first limit to overlapping portion of x,y    
    x = x(x> min(y)-10 & x < max(y)+10);
    if isempty(x)
        details.delay = 0;
        return;
    end
    y = y(y> min(x)-10 & y < max(x)+10);
    splitpt = clust.SplitTimes(x,y);
    details.ts = ts(1:end-1)+0.0001/2;
    for j = 1:size(splitpt,1)-1
        xi = x(splitpt(j,1):splitpt(j+1,1));
        yi = y(splitpt(j,2):splitpt(j+1,2));
        a = bsxfun(@minus,xi,yi');
%        xid = find(abs(a) < 0.001);
%histc(a(:)') gives the same size(1) even if length xi/yi == 1
        as{j} = histc(a(:)',ts);        
    end
    xc = sum(cat(1,as{:}),1);
    xc = xc(1:end-1);
%find max excluding end points    
    [a,b] = max(xc(2:end-1));
    b = b+1;
    details.xc = xc;
    nsync = sum(xc(b-1:b+1));
    details.wmax = max(xc([2:b-2 b+2:end])); %max count outside center bin. should be low
    xi = cat(1,as{:});
    delay = ts(b)+diff(ts(100:101))/2;
    nsyncw = sum(xc(2:end-1))-nsync; %counts over 10ms except central bin
    ratio = 3/(length(xc)-5); 
    bindur = 3/10;

elseif method == 1 %slightly quicker
    if splits > 1
        npts = floor(length(y)/splits);
        for j = 1:splits
            o = (j-1)*npts;
            if j == splits
                yi = y(o:end);
            else
                yi  = y(o + [1:npts]);
            end
            a = bsxfun(@minus,x,yi');
            xi = find(abs(a) < 0.001);
%possibly better way to estimate expected rate taht
%allows for longer timescale correlations
%            gi(j) = sum(abs(a(:)) < 0.5);
            as{j} = a(xi);
%if empty size(as{j{) is usually 0x1, and N*1 if not.             
        end
        xi = cat(1,as{:});
        delay = median(xi);
        nsync = sum(abs(xi - delay) < 0.0002);
        nsyncw = length(xi)-nsync; %counts over 1ms except central bin
    else
        a = bsxfun(@minus,x,y');
        xi = find(abs(a) < 0.001);
        delay = median(a(xi));
    end
    bindur = 2/10;
    ratio = 1/4; %duration(nsync)/duration(nsyncw)
    if splits > 1
    elseif nargout > 2
        [details.id{2}, details.id{1}] = find(abs(a - delay) < 0.0002);
        nsync = length(details.id);
        nsyncw = length(xi)-nsync;
    else
        nsync = sum(abs(a(:) - delay) < 0.0002);
        nsyncw = length(xi)-nsync;
    end
    if length(xi) == nsync %all spikes in central ms are in same bin
        if nsync == 1
            nsyncw = 1;
        else
            xc = histc(a(:)-delay,[-0.01:0.0002:0.01]);
            dx = sort(xc);
            nsyncw = dx(end-1); %second largest peak. 
        end
    end
else
    X = repmat(x(:),1,length(y));
    Y = repmat(y(:)',length(x),1);
    d = X-Y;
    nsync = sum(abs(d(:)) < 0.001);
end
details.nsync = nsync;
details.nsyncw = nsyncw;
details.delay = delay;
details.adjust = adjust;
details.rate = 1./isi;
details.wratio = ratio;
[~, details.midpt] = min(abs(ts));
ns = nsyncw * ratio; 

%details returns real counts.  But now adjust for expected rate
ansync = nsync-(adjust .* bindur);
%adjust is expect counts per ms base on overall mean rate
e = ansync./[length(x)-adjust(1)*bindur length(y)-adjust(2)*bindur];
e([3 4]) = nsyncw./[length(x)-adjust(1).*bindur/ratio  length(y)-adjust(2).*(bindur/ratio)];
if e(3) > e(1) || e(4) > e(2)% how is this poss - when adjustment for rate is main factor
%i.e. if co-incidence < expected from mean     
%    ns = nsyncw * ratio;
end
     
%adjust based on other bins in correlogram = accounts for correlations at
%longer timescales
e(7) = nsync./details.wmax; %both unadjusted
e([3 4]) = (details.wmax-ns) ./ [length(x)-ns  length(y)-ns];
e([5 6]) = (nsync-ns) ./ [length(x)-ns  length(y)-ns];
e(8) = sum(xc); %if only one co-incidence in the range spanned, can't tell if
%its real or not.  What if one bin has two, and rest are empty? etc. So Sum
%is good to know

function t = absdiff(a,b)
t = abs(a-b) < 0.001;

