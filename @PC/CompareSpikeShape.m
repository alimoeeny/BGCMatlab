function varargout = CompareSpikeShape(A,varargin)
%PC.CompareSpikeShape(C) across clusters
%C is a cell array of cluster structs
%PC.CompareSpikeShape(C, xid) across clusters
%C is a cell array of cluster structs, xid is an Nx3 matrix of [e p cl]
%triplets
%%PC.CompareSpikeShape(..., 'tag',tag) plots into named figure
%Plot MeanSpike for two 
%
%C = PC.CompareSpikeShapes(C, xid) returns a cell array of Clusters
%     can be used with PC.ShapeCorr
%
% see also PC.ShapeCorr


j = 1;
legendloc = 'LowerRight';
checkdrifts = 0;
xstr = '';
Array = [];

if ~isempty(varargin) && isnumeric(varargin{1}) && size(varargin{1},2) == 3
    xlist = varargin{1};
    for j = 1:size(xlist,1);
        C{j} = PC.GetClusterInfo(A,xlist(j,:));
    end
    try
        PC.CompareSpikeShape(C,varargin{2:end});
    catch ME
        CheckExceptions(ME);
    end
    varargout{1} = C;
    return;
end
timeshift = zeros(1,length(A));
strargs = cell2cellstr(varargin);
probes = [];
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'spacing')
        Array = varargin{j};
    elseif strncmpi(varargin{j},'bestdelay',6)
        [xc, X] = PC.ShapeCorr(A, 'delays', 5,'drifts',checkdrifts);
        timeshift(2) = X.timeshift;
    elseif strncmpi(varargin{j},'delay',3)
        j = j+1;
        timeshift  = varargin{j};
    elseif strncmpi(varargin{j},'drifts',5)
        j = j+1;
        checkdrifts = varargin{j};
    elseif strncmpi(varargin{j},'legendloc',7)
        j = j+1;
        legendloc = varargin{j};
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        GetFigure(varargin{j});
    elseif strncmpi(varargin{j},'probes',7)
        j = j+1;
        probes = varargin{j};        
    elseif strncmpi(varargin{j},'xstring',3)
        j = j+1;
        xstr = varargin{j};
    end
    j = j+1;
end


c = mycolors('ncolors',24);
linestyles = {'-' '--' ':' '-.'};
hold off;
if sum(strcmp('image',strargs));
    mysubplot(1,2,1);
    PC.PlotMeanSpike(A{1},0,A{1}.cluster,'imageonly');
    pos = get(get(gca,'title'),'extent');
    text(pos(1)+pos(3),pos(2),xstr,'VerticalAlignment','bottom');
    
    mysubplot(1,2,2);
    PC.PlotMeanSpike(A{2},0,A{2}.cluster,'imageonly');
else
    
for p = 1:length(A);
    Sa = A{p}.MeanSpike.ms;
    if isempty(probes)
        probes = 1:size(Sa,1);
    end
    for j = 1:length(probes) 
        t = [1:size(Sa,2)]-timeshift(p);
        labels{p} = sprintf('E%dP%dcl%d vs E%dP%dvs%d',A{p}.exptno,A{p}.probe,A{p}.cluster);
        h(p) = plot(t,Sa(probes(j),:),linestyles{p},'color',c{j});
        hold on;
    end
    if p > 1
        [xc(p-1), details] = PC.ShapeCorr(A{1},A{p},'delay',timeshift(p)-timeshift(1),'drifts',checkdrifts,Array);
        if details.probeshift ~= 0
            xstr = sprintf(' P+%d',details.probeshift);
        else
            xstr = '';
        end
    end
end
axis('tight');
mylegend(h,labels,legendloc);
title(sprintf('xc%.3f dt%.0f%s',xc,diff(timeshift),xstr));
end
if nargout > 0
    varargout{1} = A;
end

