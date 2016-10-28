function varargout = MakeVarMatrix(DATA, Clusters, varargin)
%PC.MakeVarMatrix(DATA, Clusters, paramlist)
% build an Nexpt X Nprobes X ncluster matrix of a cluster property
%PC.MakeVarMatrix(DATA, paramlist)
params = {};
twod = 0;
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j}, 'test')
    elseif strcmp(varargin{j}, '2D')
        twod = 1;
    elseif ischar(varargin{j})
        params = {params{:} varargin{j}};
    end
    j = j+1;
end
if nargin == 1 || isempty(Clusters)
    Clusters = getappdata(DATA.toplevel,'Clusters');
end
    if twod
        nc = 1;
    else
        nc = max(DATA.nclusters(:));
    end
    tic;
     for e = 1:length(Clusters)
         for p = 1:length(Clusters{e})
             for c = 1:nc
                 for k = 1:length(params)
                     if iscell(Clusters{e})
                         varargout{k}(e,p,c) = clust.GetValue(Clusters{e}{p},c,params{k});
                     else
                         varargout{k}(e,p,c) = clust.GetValue(Clusters{e}(p),c,params{k});
                     end
                 end
             end
         end
     end
    toc
    