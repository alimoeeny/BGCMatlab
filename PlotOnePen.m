function pen = PlotOnePen(pen,varargin)
%PlotOnePen(penlog,varargin) plots data in a penetration log named in
%penlog
%PlotOnePen(penlog,'allcomments') shows comments
%PlotOnePen(penlog,'allcomments') shows comments
%
%see also PlotComments
COMMENTS = 1;
RFS = 2;
RFS3D = 3;
plottype = COMMENTS;
plotnames = {'Comments' 'RFs' 'RF3D' 'Comments'};

holdon = 0;
if ischar(pen)
    pen = ReadPen(pen,'noplot');
end

fontsize = 18;
cmtypes = [0 1];
showcells = 1;
colors = mycolors;
j = 1;
while j <= nargin -1
    if strncmpi(varargin{j},'allcomments',6)
        j = j+1;
        cmtypes = [0 1 2];
    elseif strncmpi(varargin{j},'comments',6)
        plottype = 'commentlist';
    elseif strncmpi(varargin{j},'cmtypes',3)
        j = j+1;
        cmtypes = varargin{j};
    elseif strncmpi(varargin{j},'rfs',3)
        plottype = 2;
    elseif strncmpi(varargin{j},'hold',3)
        holdon = 1;
    elseif strncmpi(varargin{j},'rf3d',3)
        plottype = 3;
    elseif strncmpi(varargin{j},'plot',3)
        j = j+1;
        plottype = varargin{j};
    end
    j = j+1;
end

if iscellstr(pen)
    PlotOnePen(pen{1}, varargin{:});
    for j = 2:length(pen)
        PlotOnePen(pen{j},'hold', varargin{:});
    end
    return;
   
end
if ~isfield(pen,'times')
    return;
end
if holdon
    hold on;
else
    hold off;
end
F = gcf;
if isnumeric(plottype)
    plottype = plotnames{plottype};
end

if strcmp(plottype,'commentlist')
    for j = 1:length(pen.comments)
        if size(pen.cmtimestr,1) >= j
            str{j} = sprintf('%s ed%.2f %s',pen.cmtimestr(j,1:8),pen.depths(pen.cmtime(j))./1000,pen.comments{j});
        else
            str{j} = sprintf('    ed%.2f %s',pen.depths(pen.cmtime(j))./1000,pen.comments{j});
        end
        str{j} = deblank(str{j});
    end
    lst = uicontrol(F, 'Style','listbox','String', str,...
        'Callback', {@HitCommentList},...
        'Tag', 'CommentDisplay',...
        'KeyPressFcn', @KeyPressed,...
        'fontsize',fontsize,...
        'max', 2,...
        'units','norm', 'Position',[0 0.1 1 0.9]);
elseif sum(strcmp(plottype,{'PlotPen'}))
    id = find(pen.times > 0);
    t = pen.times(id);
    plot(pen.times(id),pen.depths(id));
    lh = diff(minmax(pen.depths(id)))/10;
    lw = diff(minmax(pen.times(id)))/10;
    hold on;
    for j = 1:length(pen.comments)
        id = pen.cmtime(j);
        if ismember(pen.cmtype(j),cmtypes)
            text(pen.times(id),pen.depths(id),pen.comments{j},'Rotation',90,...
                'color',colors{pen.cmtype(j)});
        end
        fprintf('%.2f(%.1f) %s',pen.depths(id),pen.times(id),pen.comments{j});
        if strncmpi(pen.comments{j},'GreyMat',5)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'k-');
        elseif strncmpi(pen.comments{j},'WhiteMat',7)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'r-');
        elseif strncmpi(pen.comments{j},'??WhiteMat',7)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'r:');
        elseif strncmpi(pen.comments{j},'?GM',3)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'k:');
        end
    end
    if isfield(pen,'enterdepth') && ~isempty(t)
    plot([t(1) t(end)],[pen.enterdepth pen.enterdepth],'k:');
    end
    if showcells
        for j = 1:length(pen.files)
            id = pen.filetimes(j);
            text(pen.times(id),pen.depths(id),pen.files{j},'Rotation',90);
        end
    end
    for j = 1:length(pen.rfs)
        t = pen.rfs(j).time;
        text(pen.times(t),pen.depths(t),sprintf('RF%.1f,%.1f',pen.rfs(j).pos(1:2)));
    end
elseif sum(strcmp(plottype,{'Comments' 'PlotPen'}))
    id = find(pen.times > 0);
    t = pen.times(id);
    plot(pen.times(id),pen.depths(id));
    lh = diff(minmax(pen.depths(id)))/10;
    lw = diff(minmax(pen.times(id)))/10;
    hold on;
    for j = 1:length(pen.comments)
        id = pen.cmtime(j);
        if ismember(pen.cmtype(j),cmtypes)
            text(pen.times(id),pen.depths(id),pen.comments{j},'Rotation',90,...
                'color',colors{pen.cmtype(j)});
        end
        fprintf('%.2f(%.1f) %s',pen.depths(id),pen.times(id),pen.comments{j});
        if strncmpi(pen.comments{j},'GreyMat',5)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'k-');
        elseif strncmpi(pen.comments{j},'WhiteMat',7)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'r-');
        elseif strncmpi(pen.comments{j},'??WhiteMat',7)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'r:');
        elseif strncmpi(pen.comments{j},'?GM',3)
            plot([pen.times(id)-lw pen.times(id)+lw],[pen.depths(id) pen.depths(id)],'k:');
        end
    end
    if isfield(pen,'enterdepth') && ~isempty(t)
    plot([t(1) t(end)],[pen.enterdepth pen.enterdepth],'k:');
    end
    if showcells
        for j = 1:length(pen.files)
            id = pen.filetimes(j);
            text(pen.times(id),pen.depths(id),pen.files{j},'Rotation',90);
        end
    end
elseif strcmp(plottype,'rfs');
    plot(0,0,'+','linewidth',2);
    hold on;
    for j = 1:length(pen.rfs)
        if isfield(pen.rfs,'depth') && ~isempty(pen.rfs(j).depth)
            z(j) = pen.rfs(j).depth;
        else
            z(j) = pen.depths(pen.rfs(j).time);
        end
        args = {};
        if isfield(pen.rfs,'fits')
            args = {args{:} 'fits',pen.rfs(j).fits};                
        end
        pos(j) = pen.rfs(j).pos(1) + i .* pen.rfs(j).pos(2);
        h(j) = plotrf(pen.rfs(j).pos,colors{j},args{:});
        if isfield(pen.rfs,'name')
            labels{j} = sprintf('%.2f %s',z(j)/1000,pen.rfs(j).name);
        else
            labels{j} = sprintf('%.2f or %.1f',z(j)/1000,pen.rfs(j).pos(5));
        end
        hold on;
    end
    lh = legend(h,labels);
    mpos = mean(pos);
    d = abs(pos-mean(pos));
    id = find(d > 1.0 .*std(d));
    if isempty(id)
        [~,id] = max(d);
    end
    if isfield(pen.rfs,'name')
        for j = 1:length(id)
            if ~isempty(pen.rfs(id(j)).name)
                text(real(pos(id(j))),imag(pos(id(j))),sprintf(pen.rfs(id(j)).name));
            end
        end
    end
    gui.PlaceLegend(lh,'TopLeft');
elseif strcmp(plottype,'RF3D');
    for j = 1:length(pen.rfs)
        if isfield(pen.rfs,'depth') && ~isempty(pen.rfs(j).depth)
            z(j) = pen.rfs(j).depth;
        else
            z(j) = pen.depths(pen.rfs(j).time);
        end
    end
    hold off;
    colors = mycolors;
    plot3(0,0,mean(z),'+');
    hold on;
    plot3(zeros(size(z)),zeros(size(z)),z,'+')
    for j = 1:length(pen.rfs)
        h = plotrf3(pen.rfs(j).pos,z(j));
        hold on;
    end
    xlabel('X');
    ylabel('Y');
    zlabel('depth');
end
title(sprintf('Pen %d at %.1f,%.1f',pen.num,pen.pos(1),pen.pos(2)));


function h = plotrf3(rf, depth)

x = [-rf(3)  rf(3) rf(3) -rf(3) -rf(3)];
y = [-rf(4)  -rf(4) rf(4) rf(4) -rf(4)];
or = rf(5) * pi/180;
xp = x .* cos(or) + y .* sin(or);
yp = -y .* cos(or) + x .* sin(or);
xp = xp+rf(1);
yp = yp+rf(2);
zp = ones(size(xp)) * depth;

h = plot3(xp,yp,zp);
arrow3([rf(1) rf(1)+rf(3) * cos(or)],[rf(2) rf(2)+rf(3) *sin(or)],zp(1:2),30,0.1,'color',get(h,'color'));

%arrow(xp,yp,or,1);

function h = plotrf(rf,color,varargin)
fits = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fits',4)
        j = j+1;
        fits = varargin{j};
    end
    j = j+1;
end

x = [-rf(3)  rf(3) rf(3) -rf(3) -rf(3)];
y = [-rf(4)  -rf(4) rf(4) rf(4) -rf(4)];
or = pi/2 -  rf(5) * pi/180;
or = rf(5) * pi/180;
xp = x .* cos(or) + y .* sin(or);
yp = -y .* cos(or) + x .* sin(or);
xp = xp+rf(1);
yp = yp+rf(2);
sz = sqrt(rf(3)^2 + rf(4)^2);

h = plot(xp,yp,'color',color);
arrow([rf(1) rf(1)+rf(3) * cos(or)],[rf(2) rf(2)+rf(3) *sin(or)],30,0.1,'color',color);
if ~isempty(fits)
    set(h,'buttondownfcn',{@HitRF, fits});
end
axis('image');

function HitRF(src, event, fits)

GetFigure('RFPlots');
hold off;
PlotExptFit(fits,'label');

function HitCommentList(a,b)
%dummy for now. Not sure what might be needed.
%currently can't use this to edit comments in pen file