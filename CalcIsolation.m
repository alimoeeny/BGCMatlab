function [x, details] = CalcIsolation(xy, idlist, clnum, varargin)
%[x, details] = CalcIsolation(pts, idlist, clnum...   calculate islolation metrics from data points
%without using GM fits
%pts is an n (points) x m (dimension of space) matrix
%idlist is a vector of n classfifications
%computes mahal separation for idlist == clnum vs idlist(~clnum)
%CalcIsolation(pts, idlist, ca, cb)
%uses idlist == ca vs idlist == cb
%x(1) is distance based on distribuion of SU points 
%x(2) is distance based on distribuion of MU points - can be negative if MU
%surrounds SU
%x(3) is like x(1) but limited to the same number of points for each
%cluster, discarding further points from cluster with most
plottype = 'none';
parentfig = 0;
ids = ismember(idlist,clnum);
id = find(ids);
nid = find(~ids);
if length(clnum) == 2
else
    nid = find(idlist ~= clnum);
end
details.overlap = NaN;
details.neardi = NaN;
x(1:4) = NaN;  %default. Be sure length matches x in good returns
plottypes = {'plothist' 'plotxy' 'plotmahal'};

j = 1;
while j <= length(varargin)
    if isnumeric(varargin{j}) && j == 1
        nid = find(ismember(idlist, varargin{j}));        
    elseif sum(strncmpi(varargin{j},plottypes,5))
        plottype = plottypes(find(strncmpi(varargin{j},plottypes,5)));
        if length(varargin) > j && isfigure(varargin{j+1})
            j = j+1;
            parentfig = varargin{j};
        end
    end
    j =j+1;
end

if length(id) < 3  || length(nid) < 3
    return;
end

%easy to get ill conditioned matrices if variables are highly correlated
wastate = warning('Off','MATLAB:illConditionedMatrix');  
wbstate = warning('Off','MATLAB:nearlySingularMatrix');  
%distance of cell points from MU dist
if length(nid) > size(xy,2) && length(id) > size(xy,2)

    %distance of SU points from MU dist
    dm = sqrt(mahal(xy(id,:),xy(nid,:)));
    
    %distance of mu points from MU dist
    dmm = sqrt(mahal(xy(nid,:),xy(nid,:)));

    %distance of MU points from SU dist
    ds = sqrt(mahal(xy(nid,:),xy(id,:)));

    %distance of SU points from SU dist
    dss = sqrt(mahal(xy(id,:),xy(id,:)));

    maxd = max([max(ds) max(dm)]);
    bins = linspace(0,maxd);
    x(1) = prctile(ds,1);
    x(1) = AllV.CalcDprime(ds,dss);
    x(2) = AllV.CalcDprime(dm,dmm);
    prctiles = prctile(dss, [50 60 70 80 90 100]);
    wts = [10 8 6 4 2 1];
    for j = 1:length(prctiles)
        c(j) = sum(ds < prctiles(j));
    end
    details.overlap = c./length(dss);
    details.overlapindex(1) = (sum(wts .* c)./sum(wts))./length(dss);
    if length(ds) > length(dss) %more mu events than su events
        [a,sid] = sort(ds); %take mu events nearest su
        details.neardi = AllV.CalcDprime(ds(sid(1:length(dss))),dss);        
        x(3) = details.neardi;
        [a,sid] = sort(dmm,'descend'); %take mu events furthrest from mu
        x(4) = AllV.CalcDprime(dm,dmm(sid(1:length(dm))));        
    else
        [a,sid] = sort(dss,'descend'); %take largeset su distances - nearest bounday
        details.neardi = AllV.CalcDprime(ds,dss(sid(1:length(ds))));        
        x(3) = details.neardi;
        [a,sid] = sort(dm);  %take smallest su events in mu dist
        x(4) = AllV.CalcDprime(dm(sid(1:length(dmm))),dmm);        
    end
else
    maxd = NaN;
end

warning(wastate);
warning(wbstate);


if ~strcmp(plottype,'none')
    if double(parentfig) > 0
        GetFigure('MahalHist','parent',parentfig);
    else
        GetFigure('MahalHist');
    end
end
if strcmp(plottype,'plothist') && ~isnan(maxd)
      
    a = hist(dss,bins);
    b = hist(ds,bins);
    c = hist(dmm,bins);
    d = hist(dm,bins);
    n(1) = max(a+b);
    n(2) = max(c+d);
    hold off;
    plot(bins,a./n(1));
    hold on;
    plot(bins,b./n(1),'g-');
    plot(bins,(a+b)./n(1),'r');
    plot(bins,c./n(2),'--');
    plot(bins,(c+d)./n(2),'r--');
    plot(bins,d./n(2),'g--');
    title(sprintf('Dprimes for SU %.2f,MU %.2f Near %.2f, %.2f',x));
elseif strcmp(plottype,'plotxy') && ~isnan(maxd)
    hold off;
    plot(xy(id,1),xy(id,2),'r.');
    hold on;
    plot(xy(nid,1),xy(nid,2),'k.');    
elseif strcmp(plottype,'plotmahal')
    hold off; 
    dm = mahal(xy,xy(nid,:));
    ds =  mahal(xy,xy(id,:));
    plot(dm(id),ds(id),'r.');
    hold on;
    plot(dm(nid),ds(nid),'k.');
end

