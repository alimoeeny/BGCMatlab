function [x, score, details] = FindEllipse(xy,idlist, varargin)
%[x, score] = FindEllipse(xy,idlist, varargin) find the ellipse that encloses points
%in xy where idlist ==2, and
%excludes points where idlist == 1
%FindEllipse(xy,idlist,'cluster',n) fits idlist == n vs idlist ~= n
% ...,'guess',X)   x is center(1) center(2) radius(1) radius(2) angle
%NB internal x(5) is in 1/10 of radians
%score(1) is # missed
%score(2) is  # false inclusions
%details.score(2) is # events that should be in
%details.score(3) = false inclusionsb

%details.score(4) is # events that should be out
guess = [];
plottype = 0;
trackfit = 0;
evalonly = 0;
ploterrors =0;
fineness = 0;
figtag = '';
j = 1;
cluster = 2;

if size(xy,1) < size(xy,2)
    xy = xy';
end

while j <= length(varargin)
    if strncmpi(varargin{j},'guess',4)
        j = j+1;
        guess = varargin{j};
    elseif strncmpi(varargin{j},'cluster',4)
        j = j+1;
        cluster = varargin{j}+1;
    elseif strncmpi(varargin{j},'fig',3)
        j = j+1;
        figtag = varargin{j};
    elseif strncmpi(varargin{j},'eval',4)
        j = j+1;
        evalonly = 1;
        guess = varargin{j};
        guess(5) = guess(5)./10;
    elseif strncmpi(varargin{j},'plot',4)
        plottype =1;
        if strncmpi(varargin{j},'ploterronly',10)
            plottype = 2;
            ploterrors = 1;
        elseif strncmpi(varargin{j},'ploterr',6)
            ploterrors = 1;
        end

    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        GetFigure(varargin{j});
    end
    j = j+1;
end
if diff(size(idlist)) > 0
    idlist = idlist';
end


%make idlist 0/1 for out/in
id = find(idlist == cluster);
nid = find(idlist ~= cluster);
idlist(id)=1;
idlist(nid) = 0;

id = id(id < size(xy,1));
if isempty(guess)
    xyr([1 2]) = mean(xy(id,:));
    xyr([3 4]) = std(xy(id,:));
    xyr(5) = 0; %angle;
else
    xyr = guess;
    xyr(5) = xyr(5)/10;
end
details.guess = xyr;
if trackfit && isappdata(0,'scores')
    rmappdata(0,'scores');
end
[~, score] = TryEllipse(xyr, xy, idlist, trackfit);
if evalonly
    x = guess;
    x(5) = x(5) .* 10;
    details.guess = [mean(xy(id,:)) std(xy(id,:)) 0];
    return;
end
options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
if trackfit
    setappdata(0,'scores',[xyr score]);
end


[x, score, details.tryangles] = TryAngles(xyr,xy, idlist, options);


quickrefine = 0;
if quickrefine
changed = 1;
while changed > 0
    [x, score, changed]  =CheckEllipseFit(x,xy,idlist, options);
end
else
[x, score, details.refine] = RefineEllipseFit(x,xy, idlist, options, fineness);
end
[~, scores] = TryEllipse(x,xy,idlist);
r = CalcClusterDistance(x,xy);
xf = x;
x(5) = x(5).*10;

details.scores = [scores(1) sum(idlist==1) scores(2) sum(idlist ==0)];
if plottype
    if plottype == 1
        hold off;
        id = find(idlist == 0);
        plot(xy(id,1),xy(id,2),'r.');
        hold on;
        id = find(idlist == 1);
        plot(xy(id,1),xy(id,2),'b.');
    end
    hold on;
    if ~isempty(guess)
        DrawEllipse(details.guess,'add','color','k');
    end
    [h, y, xx] = DrawEllipse(x,'add');
    y = repmat(y,2,1);
    y(1,:) = xx;
%    yr = CalcClusterDistance(xf,y');
%    r = clust.EllipseDistance(y',x);
    xid = find(r > 1 & idlist ==1);
    nxid = find(r < 1 & idlist == 0 );
    if ploterrors
        hold on;
        plot(xy(xid,1),xy(xid,2),'m+','markersize',6);
        plot(xy(nxid,1),xy(nxid,2),'g+','markersize',4);
    end
    title(sprintf('%d/%d missed, %d incorrectly included',length(xid),sum(idlist),length(nxid)));
    hold off; 
    
end


function [bestx, minscore, details] = TryAngles(x,xy, idlist, options)
changed = 0;
trackfit =0;
angles = [0 pi/16 pi/8 pi/4 3 * pi/8 7 * pi/16];
id = find(idlist ==1);
allx = [];
allguess = [];
scores = [];
xyr = x;
for a = angles(:)'
%In tryellipse, the data is rotated, then radius is calculated
%so use same sign for a to modiify sd guess;
    xxy = xyrotate(xy,a);
    xyr([3 4]) = std(xxy(id,:));
    xyr(5) = a/10;
    [~,b] = TryEllipse(xxy,xyr,idlist);
    [x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, xyr,options,  xy, idlist, trackfit);
    allx(end+1,:) = x;
    xyr(5) = a;
    allguess(end+1,:) = xyr;
    [~, b] = TryEllipse(x, xy, idlist, trackfit);
    scores(end+1,:) = b;
end
showguess = 0;
if showguess
    hold off;
    PlotND(xy,[],'idlist',idlist);
    hold on;
    for j = 1:size(allguess,1)
        DrawEllipse(allguess(j,:),'add');
    end
end

details.scores = scores;

[minscore,b] = min(sum(scores'));
bestx = allx(b,:);

scales = [0.6 0.8 1/0.8 1/0.6];
xyr = bestx;
for a = scales(:)'
    xyr(3) = bestx(3) * a;
    [x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, x,options,  xy, idlist, trackfit);
    [~, b] = TryEllipse(x, xy, idlist, trackfit);
    scores(end+1,:) = b;
    allx(end+1,:) = x;
end

[minscore,b] = min(sum(scores'));
bestx = allx(b,:);

x = bestx;
xyr = bestx;
for a = scales(:)'
    xyr(4) = bestx(4) * a;
    [x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, x,options,  xy, idlist, trackfit);
    [~, b] = TryEllipse(x, xy, idlist, trackfit);
    scores(end+1,:) = b;
    allx(end+1,:) = x;
end
[minscore,b] = min(sum(scores'));
bestx = allx(b,:);
details.scores = scores;
details.allx = allx;

function [bestx, minscore, details] = RefineEllipseFit(x,xy, idlist, options, fineness)
changed = 0;
[~, score] = TryEllipse(x, xy, idlist);
scores = [];
trackfit = 0;
xyr = x;
if fineness == 0
diffs = x(5) + [-0.05:0.005:+0.05];
else
diffs = x(5) + [-0.003:0.0001:+0.003];
end
bestx  = xyr;
allx = [];
for a = diffs(:)'
    xyr(5) = a;
    [x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, xyr,options,  xy, idlist, trackfit);
    allx(end+1,:) = x;
    [~, b] = TryEllipse(x, xy, idlist, trackfit);
    scores(end+1,:) = b;
end
details.scores = scores;
details.angles = diffs - bestx(5);
[minscore,b] = min(sum(scores'));
bestx = allx(b,:);
x = bestx;
x(3) = x(3) * 2;
[x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, x,options,  xy, idlist, trackfit);
if score < minscore
    bestx  = x;
end
x = bestx;
x(4) = x(4) * 2;
[x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, x,options,  xy, idlist, trackfit);
if score < minscore
    bestx  = x;
end


function [x, score, changed] = CheckEllipseFit(x,xy, idlist, options)
changed = 0;
[~, score] = TryEllipse(x, xy, idlist);
scores = [];
trackfit = 0;
diffs = -0.01:0.001:+0.01;
for a = diffs(:)'
    [~, b] = TryEllipse([x(1:4) x(5) + a], xy, idlist, trackfit);
    scores(end+1,:) = b;
end
[a,b] = min(sum(scores,2));
score = scores(b,:);
if diffs(b) ~= 0
    changed = 1;
    x(5) = x(5) + diffs(b);
    [x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, x,options,  xy, idlist, trackfit);    
    [~, score] = TryEllipse(x, xy, idlist);
else
    [x, newscore, exitflag, fitdetails] = fminsearch(@TryEllipse, x,options,  xy, idlist, trackfit);    
    if newscore < sum(score)
        changed = 1;
        [~, score] = TryEllipse(x, xy, idlist);
    end
end


function [cost, score] = TryEllipse(xyr, xy, idlist, trackfit)
%Tryellips uses radians/10 for the rotation, as this gives better
%convergence (!) 
%here idlist is true/false (1,0);
if nargin < 4
    trackfit = 0;
end

r = CalcClusterDistance(xyr,xy);
score(2) = sum(r <= 1 & ~idlist); %false inclusion
score(1) = sum(r > 1 & idlist); %missed from cluster
cost = sum(score);
if trackfit
scores = getappdata(0,'scores');
if ~isempty(scores)
    scores = [scores; xyr score];
    setappdata(0,'scores',scores);
end
end


function TestDistance(xyr,xy)
  xys = xyrotate(xy(:,1)-xyr(1),xy(:,2)-xyr(2),xyr(5));
  d = sqrt(((xys(:,1))./xyr(3)).^2 + ((xys(:,2))./xyr(4)).^2);

function r = CalcClusterDistance(xyr, xy)
  xys = xyrotate(xy(:,1)-xyr(1),xy(:,2)-xyr(2),xyr(5).*10);
  r = sqrt(((xys(:,1))./xyr(3)).^2 + ((xys(:,2))./xyr(4)).^2);
