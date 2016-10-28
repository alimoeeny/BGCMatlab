function h = myscatter(x,y,symbol,varargin)
%myscatter(x,y,symbol, ... makes a scatterplot where each datapoint has a callback
%function. Default is to print out index of data point in the vectors x,y
%
% myscatter(x,y,'.','buttonpress',@fcn)
%               sets the function called after a press
% myscatter(x,y,'.','ids',idlist)
%               sets the id numbers that are passed to @fcn
% myscatter(x,y,'.','colors',colors)
%                  sets color scheme
% myscatter(x,y,'.','colorids',colors)
%               setd color value to ues for each point. size(colors) equals size(x)
% myscatter(x,y,'.','fillids',colors)
%                sets fill color in the same way
% myscatter(x,y,'.','labels',labels)
%               attaches a label to each point in the callback
% myscatter(x,y,'.','byline') does not make separate callback for each point


fcn = @scatterhit;
idlist = [];
bidlist = [];
plotargs = {};
ptlabels = {};
colors = [];
colorids = [];
fillids = [];
callfcn = [];
byline = 0;
h = [];
if isempty(y) && size(x,2) == 2
    y = x(:,2);
    x = x(:,1);
end
if size(x,2) == 1
    x = x';
end
if size(y,2) == 1
    y = y';
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'buttonpress',10)
        j = j+1;
        fcn = varargin{j};
    elseif strncmpi(varargin{j},'byline',5)
        byline = 1;
    elseif strncmpi(varargin{j},'byid',4)
        byline = 2;
    elseif strncmpi(varargin{j},'callback',8)
        j = j+1;
        callfcn = varargin{j};
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        idlist = varargin{j};
    elseif strncmpi(varargin{j},'labels',5)
        j = j+1;
        ptlabels = varargin{j};
    elseif strncmpi(varargin{j},'colorids',8)
        j = j+1;
        colorids = varargin{j};
    elseif strncmpi(varargin{j},'fillids',7)
        j = j+1;
        fillids = varargin{j};
    elseif strncmpi(varargin{j},'colors',6)
        j = j+1;
        colors = varargin{j};
    elseif strncmpi(varargin{j},'spkid',5)
        j = j+1;
        idlist = varargin{j};
        byline = 1;
        colors = mycolors('spkcolors');
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end

if ~isempty(colorids) && isempty(colors)
    colors = mycolors;
end
if isempty(idlist)
    idlist = ones(size(x));
elseif size(idlist,2) == 1
    idlist = idlist';
end
if isempty(bidlist)
   bidlist = 1:size(x,1);
end
np = 0;

%put the tuning curve data that generated X,Y into DATA
DATA.x = x;
DATA.y = y;
if byline == 1
    np = 1;
    for k = 1:size(x,1)
        ln = unique(idlist);
        for j = 1:length(ln)
            id = find(idlist == ln(j));            
            h(np) =  plot(x(k,id),y(k,id),symbol,'buttondownfcn',{fcn, ln(j), bidlist(k)},plotargs{:});
            hold on;
            if ~isempty(colorids)
                set(h(np),'color',colors{colorids(k,j)});
            elseif ~isempty(colors)
                set(h(np),'color',colors{k,ln(j)});
            else
        end
        np = np+1;
        end
    end
elseif byline == 2
    ln = unique(idlist);
    np = 1;
    hold off;
    for j = 1:length(ln)
        id = find(idlist == ln(j));
        h(np) =  plot(x(id),y(id),symbol,'buttondownfcn',{fcn, ln(j), ln(j)},plotargs{:});
        hold on;
        np = np+1;
    end    
else
    if isempty(fillids)
        fillids = zeros(size(x));
    end
    plotpt = ones(size(x));
    plotpt(colorids==0) = 0;
    for ia = 1:size(x,3)
    for k = 1:size(x,1)
        for j = 1:size(x,2)
            if plotpt(k,j,ia)
            np = np+1;
            
            %can also plot all teh data in one line, but then the callback function has
            %to do work to find the data point closest to the mouse.  For large
            %datasets, its faster to plot just one line.
            h(np) =  plot(x(k,j,ia),y(k,j,ia),symbol,'buttondownfcn',{fcn, idlist(:,j), bidlist(k)},plotargs{:});
            if ~isempty(colorids)
                set(h(np),'color',colors{colorids(k,j,ia)});
            elseif ~isempty(colors)
                set(h(np),'color',colors{k,j,ia});
            end
            if fillids(k,j,ia)
                set(h(np),'markerfacecolor',colors{fillids(k,j,ia)});
            end
            hold on;
            end
        end
    end
      ln = ia;
    end
end
setappdata(gcf,'XYData',DATA);
setappdata(gcf,'callfcn',callfcn);
if ~isempty(ptlabels)
    setappdata(gcf,'Labels',ptlabels);
end


function scatterhit(a,b,ida, idb)

pos = get(gca,'currentpoint');
labels = getappdata(gcf,'Labels');
DATA = getappdata(gcf,'XYData');
if isempty(labels)
    lbl = '';
elseif sum(size(labels) > 1)  == 2
    lbl = labels{idb,ida};    
elseif ida == 1
    lbl = labels{idb};
else
    lbl = labels{ida};
end
%If there was app data stoerd in the figure, would be able to  
%use that to E.g. plot tuning curves. 
fprintf('Point %s, %d at %s %s %s\n',sprintf('%d,',ida), idb,num2str(pos(1,1)),num2str(pos(1,2)),lbl);
callfcn = getappdata(gcf,'callfcn');
if ~isempty(callfcn)
    if iscell(callfcn)
        feval(callfcn{1}, a,b,idb,ida,callfcn{2:end},'figure','Tagb');
    else
        feval(callfcn, a,b,idb,ida,'figure','Tagb');
    end
end
