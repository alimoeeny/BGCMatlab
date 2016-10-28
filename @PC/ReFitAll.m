function DATA = ReFitAll(DATA, fittype)    Clusters = getappdata(DATA.toplevel,'Clusters');    Expts = getappdata(DATA.toplevel,'Expts');    if strcmp(fittype,'sdindex')        DATA.xysdindex = [];    end    ts = now;    for j = 1:length(Clusters)        DATA.Expt = Expts{j};        for k = 1:length(Clusters{j})            C = Clusters{j}{k};            if strcmp(fittype,'3means')             a = GMfit(C.xy, 3,1);            DATA.mahal3(j,k) = GMdprime(a);            elseif strcmp(fittype,'sdindex')                sdx = PC.PlotXYSequence(DATA, [j k],'noplot');                DATA.xysdindex(j,k) = sdx(end);            end        end        fprintf('Ex%d %.2f\n',j,mytoc(ts));    end    PC.SaveExtras(DATA);    set(DATA.toplevel,'UserData',DATA);