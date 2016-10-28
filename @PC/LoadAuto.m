function [DATA,A] = LoadAuto(DATA, varargin)
%DATA = PC.LoadAuto(DATA) Loads auto cluster files
%DATA = PC.LoadAuto(DATA,'ecker') Loads files from ecker subdir
%DATA = PC.LoadAuto(DATA,'ecker','setdata') sets the 'AutoClusters' appdata
%                         with the loaded Clusters (merged with any existing

DATA = GetDataFromFig(DATA);
setdata = 0;
autotype = 'gmm';
needdetails = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'ecker',5)
        autotype = 'ecker';
    elseif strncmpi(varargin{j},'setdata',5)
        setdata = 1;
    end
    j = j+1;
end


args = {};
if ischar(DATA) %manuual loading
    ddir = DATA;
    clear DATA;
    DATA.datadir = ddir;
end

if strcmp(autotype,'ecker')
    ddir = [DATA.datadir '/Ecker'];
else
    ddir = DATA.datadir;
end
DATA.autodatadir = ddir;

if ~isfield(DATA,'ArrayConfig')
    DATA.ArrayConfig = GetArrayConfig(DATA.datadir);
end

d = mydir([ddir '/Expt*AutoClusterTimes.mat']);
ts = now;

if needdetails
    args = {args{:} 'needdetails'};
end
for j = 1:length(d)
    [C{j},F{j}, D{j} CD{j}] = LoadCluster(d(j).name,'getxy','autoonly','savespace',args{:});
    expts(j) = GetExptNumber(C{j});
    if isfield(D{j},'loaddur')
        loadsize(j) = sum(D{j}.loadbytes)./(1024 * 1024);
        loadrate(j) = loadsize(j)./sum(D{j}.loaddur);
        if loadrate(j) < 20
            fprintf('Load Rate for %s only %.1f MB/sec\n',d(j).name,loadrate(j));            
        end
    else
        loadrate(j) = NaN;
    end
end
if isempty(d)
    return;
end
fprintf('LoadAuto: Mean load rate %.2f. %.1fMb Real time %.1f\n',nanmean(loadrate),sum(loadsize),mytoc(ts));
%now sort by expt nnumber
[expts,b] = sort(expts);
C = C(b);
D = D(b);
CD = CD(b);

[C, GM] = PC.CondenseClusters(C);
A = getappdata(DATA.toplevel,'AutoClusters');
if isempty(A)
    A = C;    
else
    for j = 1:length(d)
        aexpts(j) = GetExptNumber(A{j});
    end
    for j = 1:length(expts)
        id = find(aexpts == expts(j));
        if length(id) == 1
            A{id} = C{j};
        elseif isempty(id)
            fprintf('Adding Autoclusters for New Expt %.1f\n',expts(j));
            A{end+1} = C{j};
        else
            fprintf('Expt %.1f is repeated in AutoClusters\n',expts(j));
        end
    end
end
for j = 1:length(A)
    ts = CellToMat(A{j},'savetime');
    fprintf('Expt%d Cluster Saves %s to %s\n',j,datestr(min(ts(:))),datestr(max(ts(:))));
    for k = 1:length(A{j})
        A{j}{k}.exptid = j;
        CD{j}{k}.exptid = j;
        if ~isfield(A{j}{k},'sign')
            A{j}{k}.sign = 0;
        end
    end
end
if setdata
    setappdata(DATA.toplevel,'AutoClusters',A);
    setappdata(DATA.toplevel,'AutoClusterInfo',D);
    setappdata(DATA.toplevel,'AutoClusterFits',GM);
    if needdetails
        setappdata(DATA.toplevel,'AutoClusterDetails',CD);
    end
    if nargout == 0
        SetData(DATA);
    end
end
    
