function ButtonPressed(src, data)DATA = GetDataFromFig(src);DATA.ts = now;start = get(gca,'CurrentPoint');DATA.elmousept.newcluster = 0;if AllV.InGraph(start,gca)    get(gcf,'selectiontype');    mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});    if DATA.profiling        DATA.elmousept.downtime = now;        fprintf('Mode is %s\n',get(gcf,'SelectionType'));    end%    DATA.elmousept.mode = mode;%only check for in ellipse if in correct space    C = AllV.ClusterInfo(DATA);    DATA = AllV.SetSpaceFromAxis(DATA, gca);    axcl = getappdata(gca,'axiscluster');    distances = [];    for j = 1:length(axcl)        CC = PC.GetClusterInfo(DATA.cluster,axcl(j));        distance = AllV.DistanceToEllipse(CC,start(1,1:2));        edgedistances(j) = abs(distance(1)-1);        distances(j) = distance(1);        if DATA.profiling            fprintf('Cl %d D%.2f,%.2f\n',axcl(j),distance);        end    end    [a,b] = min(distances);    if a < 1  && axcl(b) ~= DATA.elmousept.cluster%force selection of this cluster        DATA.elmousept.newcluster = axcl(b);        DATA.elmousept.shape = 0;        SetData(DATA);        if DATA.profiling            fprintf('Inside Cl  %d\n',axcl(b));        end    end    DATA.elmousept.start = start(1,1:2);    distance = AllV.DistanceToEllipse(DATA.elmousept,start(1,1:2));%Not so simple. When an ellipse is drawn in a new space, C.space is not%updated until clicking inside. %    if length(DATA.elmousept.pcplot) == 2 && length(C.space) > 1%        if isfield(C,'space') &&  sum(C.space(end-1:end) == DATA.elmousept.pcplot) < 2%            distance = 100%        end    %end        if ishandle(DATA.elmousept.h)        if get(DATA.elmousept.h,'Parent') == gca  %only inside if its teh right subplot            axisok = 1;        else            axisok = 0;        end    else        axisok = 1;    end    if ismember(DATA.elmousept.cluster,axcl) %this cluster ellipse is in this axis        axisok = 1;    else        axisok = 1;    end    if DATA.profiling && isfield(DATA.elmousept,'markpt');        fprintf('%.1f,%.1f Distance %.1f,%.1f mark at %.1f,%.1f\n',start(1,1:2),distance,DATA.elmousept.markpt);    end    yl = get(gca,'ylim');    xl = get(gca,'xlim');    set(gca,'xlimmode','manual','ylimmode','manual'); %dont rescale for ellispes    DATA.elmousept.aspectratio = diff(yl)./diff(xl);    if  mode  == 2 % R button press. Inside an existing ellispe sets cutting to that space        DATA.elmousept.down = 3;        DATA.elmousept.done = 0;        [x, c] =  AllV.FindNearestCluster(DATA, start(1,1:2));        if x < 1            DATA = AllV.SetEllipseDrawing(DATA,0,'cluster',c);        else            mode = 3;            DATA.elmousept.down = 5;            DATA.elmousept.pts = start(1,1:2);        end    elseif mode == 3 %shift press        fprintf('Shift is Down DIstances %.2f,%.2f\n',distance);        DATA.elmousept.down = 5;        DATA.elmousept.pts = start(1,1:2);    elseif distance(2) < 0.2        DATA.elmousept.down = 4;        fprintf('Setting Angle\n');    elseif distance(1) < 1.05 && axisok%test fails for NaN%pressed inside ellipse, so just tmove it.         DATA.elmousept.down = 2;        DATA.elmousept.pos =[0  0 0 0 ];        DATA.elmousept.start = start(1,1:2);        DATA.elmousept.axis = gca;        %set up pos in case released with no drag        DATA.elmousept.pos(1) = DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);        DATA.elmousept.pos(2) = DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);        DATA.elmousept.pos(3) = DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);        DATA.elmousept.pos(4) = DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);    else        DATA.elmousept.down = 1;        DATA.elmousept.done = 0;        DATA.elmousept.pos =[start(1,1) start(1,2) 0 0 ];        DATA.elmousept.axis = gca;        yl = get(gca,'ylim');        xl = get(gca,'xlim');        DATA.elmousept.aspectratio = diff(yl)./diff(xl);  %if drawing in a new space, set cluster angle to zero.            DATA.elmousept.pcplot = get(gca,'UserData');        if length(DATA.elmousept.pcplot) == 2 && length(C.space) > 1            if isfield(C,'space') &&  sum(C.space(end-1:end) == DATA.elmousept.pcplot) < 2                DATA.elmousept.angle = 0;            end        end                if 0        if ~isfield(DATA.elmousept,'angle')            DATA.elmouept.angle = 0 ;        end        if ~isfield(DATA.elmousept,'plotargs')            DATA.elmouept.plotargs = {} ;        end        end    endset(DATA.toplevel,'UserData',DATA);end