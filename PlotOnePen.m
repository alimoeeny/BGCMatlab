function PlotOnePen(pen,varargin)

COMMENTS = 1;
RFS = 2;
RFS3D = 3;
plottype = COMMENTS;

if ischar(pen)
    pen = ReadPen(pen,'noplot');
end

cmtypes = [0 1];
showcells = 1;
colors = mycolors;
j = 1;
while j <= nargin -1
    if strncmpi(varargin{j},'cmtypes',3)
        j = j+1;
        cmtypes = varargin{j};
    elseif strncmpi(varargin{j},'plot',3)
        j = j+1;
        plottype = varargin{j};
    end
    j = j+1;
end

if ismember(plottype,[COMMENTS])
    plot(pen.times,pen.depths);
    lh = range(pen.depths)/10;
    lw = range(pen.times)/10;
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
    if isfield(pen,'enterdepth')
    plot([pen.times(1) pen.times(end)],[pen.enterdepth pen.enterdepth],'k:');
    end
    if showcells
        for j = 1:length(pen.files)
            id = pen.filetimes(j);
            text(pen.times(id),pen.depths(id),pen.files{j},'Rotation',90);
        end
    end
elseif ismember(plottype,[RFS]);
    hold off;
    plot(0,0,'+','linewidth',2);
    hold on;
    for j = 1:length(pen.rfs)
        z(j) = pen.depths(pen.rfs(j).time);
        h(j) = plotrf(pen.rfs(j).pos,colors{j});
        labels{j} = sprintf('%.2f or %.1f',z(j)/1000,pen.rfs(j).pos(5));
        hold on;
    end
    legend(h,labels);
elseif ismember(plottype,[RFS3D]);
    for j = 1:length(pen.rfs)
        z(j) = pen.depths(pen.rfs(j).time);
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
end
title(sprintf('Pen %d',pen.num));


function h = plotrf3(rf, depth)

x = [-rf(3)  rf(3) rf(3) -rf(3) -rf(3)];
y = [-rf(4)  -rf(4) rf(4) rf(4) -rf(4)];
or = rf(5) * pi/180;
xp = x .* cos(or) + y .* sin(or);
yp = y .* cos(or) - x .* sin(or);
xp = xp+rf(1);
yp = yp+rf(2);
zp = ones(size(xp)) * depth;

h = plot3(xp,yp,zp);
%arrow(xp,yp,or,1);

function h = plotrf(rf,color)

x = [-rf(3)  rf(3) rf(3) -rf(3) -rf(3)];
y = [-rf(4)  -rf(4) rf(4) rf(4) -rf(4)];
or = pi/2 -  rf(5) * pi/180;
xp = x .* cos(or) + y .* sin(or);
yp = -y .* cos(or) + x .* sin(or);
xp = xp+rf(1);
yp = yp+rf(2);
sz = sqrt(rf(3)^2 + rf(4)^2);

h = plot(xp,yp,'color',color);
arrow([rf(1) rf(1)+2 * cos(or)],[rf(2) rf(2)+2 *sin(or)],30,0.3,'color',color);
axis('image');