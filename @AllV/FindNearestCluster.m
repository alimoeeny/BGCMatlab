function [distance, cluster] = FindNearestCluster(DATA, pos)    distance = NaN;    cluster = 0;        pcplot = DATA.elmousept.pcplot;    if length(pcplot) == 3        pcplot = pcplot(2:3);    end    if DATA.cluster.space(1) == DATA.elmousept.pcspace && sum(pcplot == DATA.cluster.space(2:end)) == length(pcplot)        d(1) = AllV.DistanceToEllipse(DATA.cluster, pos,'radius');    else        d(1) = NaN;    end    for j = 1:length(DATA.cluster.next)        if ~isempty(DATA.cluster.next{j}) && DATA.cluster.next{j}.space(1) == DATA.elmousept.pcspace && sum(pcplot == DATA.cluster.next{j}.space(2:end)) == length(pcplot)            d(j+1) = AllV.DistanceToEllipse(DATA.cluster.next{j}, pos,'radius');        else            d(j+1) = NaN;        end    end    [distance, cluster] = min(d);    if DATA.profiling        for j = 1:length(d)            fprintf('D%d %.1f\n',d(j));        end    end    