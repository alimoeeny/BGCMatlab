function need = ClusterNeedsRefresh(C, I)        need = 0;    if ~isfield(I,'loadname') %missing info, can't tell        nee = NaN;        return;    end    d = dir(I.loadname);    if strfind(I.loadname,'AutoCluster')        name = strrep(I.loadname,'AutoCluster','Cluster');        if exist(name)            need = 3;        end    end    if d.datenum > I.loadtime        need = 4;    end    if need > 0        d = dir(strrep(I.loadname,'ClusterTimes','ClusterTimesDetails'));        if isempty(d) || d.datenum < I.loadtime            need = 2; %Details not modified = quick save.          end    end