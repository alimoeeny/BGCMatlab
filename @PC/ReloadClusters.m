 function Clusters = ReloadClusters(DATA, eid)      cfile = [DATA.name '/' DATA.strings{eid}]; if strfind(cfile,'AutoClusterTimes')     Clusters = LoadCluster(cfile,'getxy'); else     afile = strrep(cfile,'ClusterTimes','AutoClusterTimes')     AutoClusters = LoadCluster(afile,'getxy');     Clusters = LoadCluster(cfile,'getxy');     for j = 1:length(AutoClusters)         if j > length(Clusters) || ~isfield(Clusters{j},'mahal') || ...                 (Clusters{j}.auto == 1 && AutoClusters{j}.savetime(1) > Clusters{j}.savetime(1))             Clusters{j} = AutoClusters{j};         end     end end for j = 1:length(Clusters)     Clusters{j}.exptid = eid; end Clusters = PC.FixClusters(Clusters);     