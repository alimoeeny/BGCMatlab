function p = GetClusterNumber(DATA, varargin)
%p = GetClusteNumber(DATA) extracts cluster # from struct/string
%if DATA is a cell array, p is a vector of cluster ids
p = NaN; 
if ischar(DATA);
    p = GetClusterFromName(DATA, varargin{:});
elseif iscell(DATA)
    for j = 1:length(DATA)
        p(j) = GetClusterNumber(DATA{j});
    end
elseif isfield(DATA,'trueprobe') && DATA.trueprobe > 0
    p = DATA.trueprobe;
elseif isfield(DATA,'probelist')
    if DATA.probe(1) > length(DATA.probelist)
        p = DATA.probe(1);
    elseif DATA.probe(1) == 0
        p = DATA.probelist(1);
    else
        p = DATA.probelist(DATA.probe(1));
    end
elseif isfield(DATA,'probe')
    p = DATA.probe(1);
elseif isfield(DATA,'Header')
    H = DATA.Header;
    if isfield(H,'clusterid')
        p = H.clusterid;
    elseif isfield(H,'Clusters')
        p = H.Clusters{1}.cluster'
    end
else
    p = NaN;
end
    

function p = GetClusterFromName(name)

p = NaN; %for now