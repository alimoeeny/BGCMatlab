function [latency, details] = sdflatency(sdf, presamples, varargin)

%[latency, details] = sdflatency(sdf, presamples, varargin) estmate latency
% of sdf. 
% returns a latency estimate from a spike density function, sdf
% presamples is the number of samples in the sdf guaranteed not to
% contain a response.
%
%sdflatency(sdf, presamples,'offset',offset)
%uses the range offset:offset+pressamples as the preperiod
%
%sdflatency(sdf, presamples,'baselie',b)


baseline = NaN;
offset = 0;
setmaxtime = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'baseline',4)
        j = j+1;
        baseline = varargin{j};
    elseif strncmpi(varargin{j},'maxt',4)
        j = j+1;
        setmaxtime = varargin{j};
    elseif strncmpi(varargin{j},'offset',4)
        j = j+1;
        offset = varargin{j};
    end
    j = j+1;
end

details.fitted = 1./zeros(size(sdf)); %NaNs

if presamples+offset > length(sdf)
    latency = NaN;
    return;
end

if size(sdf,2) == 1
    sdf = sdf';
end

t(1) = presamples+offset;
prerate = mean(sdf(1+offset:presamples+offset));
if sum(~isnan(sdf)) == 0
    latency = NaN;
    return;
end
if setmaxtime
    [maxv, maxt] = max(sdf(1:setmaxtime));
else
    [maxv, maxt] = max(sdf);
end

if maxt < presamples
    details.premax = maxv;
    [maxv, maxt] = max(sdf(presamples:end));
    maxt = maxt+presamples-1;
end
details.maxrate = maxv;
details.tmax = maxt;
halfv = prerate + (maxv-prerate).*0.8;
id = find(sdf(1:maxt) > halfv);
if length(id) == 1 && id(1) == maxt
    pt = maxt;
else
    pt = id(2); %first point that reaches 80% max
    if pt < presamples
        id = find(sdf(presamples:maxt) > halfv);
        if length(id) > 1
            id = id + presamples-1;
            pt = id(1);
        else
            latency = NaN;
            details.failure = 'MaxV is only point > 80%';
            details.baseline = 0;
            return;
        end
    end
end
halfv = prerate + (maxv-prerate)/2;
id = find(sdf(1:pt) < halfv);
if isempty(id) %never less than halfmax. Means exterything before max
end
halft = id(end);

%This is a key choice. For strong sdfs with a curved onset/ noise baseline
%it pays to use a hihger value for v.  For irregular repsosnes, setting a 
%high value for v can give very late latencies
v = prerate + (maxv-prerate)/4;



id = find(sdf(1:halft) < v);
if isempty(id)
    vt = halft;
else
    if maxt-id(end) < 20 && maxv > std(sdf) .* 10 
        v = prerate + (maxv-prerate)/2;
        id = find(sdf(1:maxt) < v);
    end        
    vt = id(end);
end

%guess latency
if vt < halft
    x(1) = vt - (halft-vt)/1.5;
    x(2) = (halfv-v)./(halft-vt); %guess slope
else
    dt = round((maxt-halft)/2);
    if dt > halft/2
        dt = halft/2;
    end
    x(2) = (diff(sdf(halft+[-dt dt])))/(2.*dt); %guess slope
    x(1) = round(halft - (sdf(halft)-prerate)./x(2));
    vt = vt-1;
end
if x(1) < t(1)
    x(1) = t(1);
end

if isnan(baseline)
    x(3) = prerate;
else
    state.baseline = baseline;
end
if halft < t(1)+presamples/10;
    halft = t(1)+presamples/2;
end

state.sdf = sdf(1:halft);
options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
[fittedparams, fval,exitflag, output] = fminsearch(@FitStart, x,options, state);
[a, b] = FitStart(fittedparams,state);
if ~isempty(b)
    details.fitted = b.fit;
    details.baseline = b.base;
else
    details.baseline = NaN;
end
latency = fittedparams(1);



function [ssq, details] = FitStart(x, state)

details = [];
x(1) = round(x(1));
if x(1) >= length(state.sdf) || x(1) < 1
    ssq = NaN;
    return;
end
if length(x) == 3
    base = x(3);
else
    base = state.baseline;
end
pred(1:x(1)) = base;
pred(x(1):length(state.sdf)) = base + (x(2) .* [0:length(state.sdf)-x(1)]);
diffs = state.sdf-pred;
ssq = sum(diffs.^2);
details.fit = pred;
details.base = base;
