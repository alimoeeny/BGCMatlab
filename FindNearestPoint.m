function [exy, id, distance] = FindNearestPoint(h, xy, varargin)
%[exy, id, distance] FindNearestPoint(h, xy, varargin) find nearest data point in line to
%current pos given by xy id returns the index in xy
%[exy, id] = FindNearestPoint(h) uses get(gca,'currentpos') for xy
%[exy, id] = FindNearestPoint() Checks all data lines in gca


distance = NaN;
if nargin == 0 || strcmp(class(h),'matlab.graphics.axis.Axes')
    if nargin == 0
        ax = gca;
    else
        ax = h;
    end
    h = findobj(allchild(ax),'flat','Type','Line');
    xy = get(gca,'CurrentPoint');
    for j = 1:length(h)
        [exy(j,:), id(j),d(j)] = FindNearestPoint(h(j), xy);
    end
    [distance,j] = min(d);
    exy = exy(j,:);
    id = id(j);
   return;     
end
if nargin ==1
    xy = get(gca,'CurrentPoint');
end
xx = xy(1);
yy = xy(1,2);
x = get(h,'xdata');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
if strcmp(get(gca,'xscale'),'log')
    xl = log(xl);
    x = (log(x)-xl(1)) ./ diff(xl);
    xx = log(xx);
else
    x = (x-xl(1))./diff(xl);
end
y = get(h,'ydata');
y = (y-yl(1))./diff(yl);
d = ((xx-xl(1))./diff(xl) + i*(yy -yl(1))./diff(yl)) -(x + i .* y);
[dd, id] = min(abs(d));
exy(1) = (x(id).*diff(xl))+xl(1);
exy(2) = (y(id).*diff(yl))+yl(1);
if strcmp(get(gca,'xscale'),'log')
    exy(1) = exp(exy(1));
end