 function PCButtonDragged(src, data)DATA = GetDataFromFig(src);x = [];if isfield(DATA,'elmousept')    if DATA.profiling && DATA.elmousept.down > 0        fprintf('Mode %d\n',DATA.elmousept.down);    endif  DATA.elmousept.down == 1    start = get(gca,'CurrentPoint');    DATA.elmousept.pos(3) = start(1,1);    DATA.elmousept.pos(4) = start(1,2);    [DATA.elmousept.h, x] = AllV.DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});%    drawnow;%    mytoc(DATA.ts);elseif  DATA.elmousept.down == 3 %set radius with R mouse button    start = get(gca,'CurrentPoint');    r = abs(start(1,1:2) - DATA.elmousept.xyr(1:2));    sina = sin(DATA.elmousept.angle);    cosa = cos(DATA.elmousept.angle);    r = r * [cosa sina; -sina cosa];    DATA.elmousept.xyr([3 4]) = r;    DATA.elmousept.pos(1) = DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);    DATA.elmousept.pos(2) = DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);    DATA.elmousept.pos(3) = DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);    DATA.elmousept.pos(4) = DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);    [DATA.elmousept.h, x] = AllV.DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});elseif  DATA.elmousept.down == 2 %moving ellipse    start = get(gca,'CurrentPoint');    DATA.elmousept.pos(1) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);    DATA.elmousept.pos(2) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);    DATA.elmousept.pos(3) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);    DATA.elmousept.pos(4) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);    [DATA.elmousept.h, x] = AllV.DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});elseif  DATA.elmousept.down == 4 %set angle    start = get(gca,'CurrentPoint');    angle = atan2(start(1,1)-DATA.elmousept.xyr(1),start(1,2)-DATA.elmousept.xyr(2));    DATA.elmousept.angle = angle;    [DATA.elmousept.h,x] = AllV.DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});elseif  DATA.elmousept.down == 5 %Draw Arbitray shape    start = get(gca,'CurrentPoint');    DATA.elmousept.pts(end+1,:) = start(1,1:2);    if ishandle(DATA.elmousept.h) && strcmp(get(DATA.elmousept.h,'type'),'line') && ...            get(DATA.elmousept.h,'parent') == gca        set(DATA.elmousept.h,'xdata',DATA.elmousept.pts(:,1));        set(DATA.elmousept.h,'ydata',DATA.elmousept.pts(:,2));    else        DATA.elmousept.h = line(DATA.elmousept.pts(:,1),DATA.elmousept.pts(:,2));    endendDATA.elmousept = CopyFields(DATA.elmousept,x,'markpt');set(DATA.toplevel,'UserData',DATA);end