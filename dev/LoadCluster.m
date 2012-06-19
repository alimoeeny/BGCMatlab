function [AllClusters, AllFullVData] = LoadCluster(dirname, eid, varargin)
%LoadCluster(dirname, expts, ...)
%Load Cluster info for an expt or list of expts
%if length(expts) = 1 returns a 1 x nprobes cell matrix, else 
%a 1 x length(expts) matrix each containing 1 x nprobes cells
%
%LoadCluster Combines info from Clusters and ClusterDetails files
%so that clst, t and Evec fields are in clusters
%LoadCluster(dirname, expts, 'gextxy') also inlcudes xy
%LoadCluster(dirname, expts, 'rawxy') also inlcudes xy with any roation
%removed (PlotClusters needs it this way
%LoadCluster(dirname, expts, 'alltimes') 
%replaces Clusters{}.times with ClusterDetails{}.t, saving memory
%(Also needed by PlotClusters)

AllClusters = {};
AllFullVData = {};
f = {'Evec' 'clst' 't'}; %fields to copy
rawxy = 0;
alltimes = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'getxy',5)
        f = {f{:} 'xy'};
    elseif strncmpi(varargin{j},'rawxy',5)
        rawxy = 1;
    elseif strncmpi(varargin{j},'alltimes',5)
        alltimes = 1;
    end
    j = j+1;
end

for j = 1:length(eid)
    name = [dirname '/Expt' num2str(eid(j)) 'ClusterTimes.mat'];
    dname = [dirname '/Expt' num2str(eid(j)) 'ClusterTimesDetails.mat'];
    daname = [dirname '/Expt' num2str(eid(j)) 'AutoClusterTimesDetails.mat'];
    aname = [dirname '/Expt' num2str(eid(j)) 'AutoClusterTimes.mat'];
    if exist(aname,'file')
        load(aname);
        AutoClusters = Clusters;
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
    else 
        AutoClusters = {};
    end
    if exist(name,'file')
        load(name);
        AllClusters{j} = Clusters;
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
        for k = 1:length(AutoClusters)
            if k > length(Clusters) || isempty(Clusters{k})
                AllClusters{j}{k} = AutoClusters{k};
            end
        end
    elseif ~isempty(AutoClusters)
        AllClusters{j} = AutoClusters;
        fprintf('Can''t read %s\n',name);
    else
        fprintf('Can''t read %s or %s\n',name,aname);
    end
    if exist(daname,'file')
        load(daname);
        AutoClusterDetails = ClusterDetails;
    else
        AutoClusterDetails = [];
    end
    if exist(dname,'file')
        load(dname);
        for k = 1:length(AutoClusterDetails)
            if k > length(ClusterDetails) || isempty(ClusterDetails{k})
                ClusterDetails{k} = AutoClusterDetails{k};
            end
        end
    else
        ClusterDetails = AutoClusterDetails;
    end
    for k = 1:length(ClusterDetails)
        for n = 1:length(f)
            if isfield(ClusterDetails{k},f{n})
                AllClusters{j}{k}.(f{n}) = ClusterDetails{k}.(f{n});
            end
        end
        if rawxy
            xy = ClusterDetails{k}.xy;
            C = AllClusters{j}{k};
            if C.shape == 0
                AllClusters{j}{k}.xy = xy;
            else
                AllClusters{j}{k}.xy = xyrotate(xy(:,1),xy(:,2),-C.angle);
            end
        end
        if alltimes && isfield(AllClusters{j}{k},'t')
            AllClusters{j}{k}.times = AllClusters{j}{k}.t;
            AllClusters{j}{k} = rmfield(AllClusters{j}{k},'t');
        end
    end
end

if length(AllClusters) == 1
    AllClusters = AllClusters{1};
    if length(AllFullVData)
    AllFullVData = AllFullVData{1};
    end
end