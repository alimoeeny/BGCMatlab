function s = str(X, varargin)
%s = str(X, varargin) make string identifying X
%typically reports Expt/Probe number for structures where this is
%valid
%also works on ArrayConfig
%str(X, 'nocluster') omits deiscription of cluster # 
%str(X, 'all') if X is a cellarray, returns an array of strings
strargs = {};
showcluster = 0;
nocluster = 0;
j = 1;
while j <= length(varargin)
    if ischar(varargin{j})
        strargs = {strargs{:} varargin{j}};
    end
j = j+1;
end
if sum(strncmp('nocluster',strargs,4))
    nocluster = 1;
elseif sum(strncmp('double',strargs,4))
end

if isnumeric(X)
    if length(X) ==1
        s = sprintf('%d',X);
    elseif length(X) == 2
        s = sprintf('%.1f',floor(X(1))+X(2)/10);
    elseif length(X) == 3
        s = sprintf('%.2f',floor(X(1))+X(2)/10+X(3)/100);
    end
elseif isfield(X,'progname') && strncmp(X.progname, 'AllVPCs',7)
    s = AllV.IDStr(X, varargin{:});
    return;
elseif isfield(X,'spacing') && isfield(X,'src') %An ArrayConfig
    s = '';
    if isfield(X,'type')
        if strcmp(X.type,'unknown') && isnan(X.spacing)
            s = 'unknown array';
        else
            s = [s X.type ' ' num2str(X.spacing) 'u'];
        end
    else
        s = [num2str(X.spacing) 'u'];
    end    
elseif iscell(X)
    alls = '';
    for j = 1:length(X)
        s{j} = str(X{j}); 
        alls = [alls s{j} ' '];
    end
    if ~sum(strcmp('all',strargs))
        s = alls;
    end
    return;
else
    showcluster = 0;
    p = GetProbeNumber(X);
    e = GetExptNumber(X);
    if isfield(X,'cluster') && X.cluster > 1  || iscluster(X)
        if nocluster == 0
            showcluster = 1;
        end
    end
    if isnan(e)
        s = sprintf('E?P%d',p);
    elseif length(p) >2
        s = sprintf('E%dP%d-%d',e,minmax(p));
    else
        s = sprintf('E%dP%d',e,p);
    end
end

if showcluster || sum(strncmp('cluster',strargs,2))
    if isfield(X,'currentcluster')
        s = sprintf('%scl%d',s,X.currentcluster)
    elseif isfield(X,'cluster') && isnumeric(X.cluster)
        s = sprintf('%s.%d',s,X.cluster);
    end
end