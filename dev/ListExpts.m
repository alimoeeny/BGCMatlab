function [names, details] = ListExpts(E, varargin)
%[names, details] = ListExpts(E, ...)   lists names of Expts in  E
% E can be a cell array of expts, or a directory name
if iscellstr(E) %list of names
    FindExpts(E, varargin);
elseif iscell(E)
    if iscellstr(E{1})
        [names, details] = FindExpts(E, varargin{:});
    else
        [names, details] = ListExptCells(E, varargin{:});
    end
elseif isdir(E)
    [names, details] = ListExptDir(E, varargin);
elseif ischar(E)
    [names, details] = ListExptDirs(E, varargin);
end

function [names, ids]= FindExpts(namelist, varargin)
names = {};

if isempty(namelist)
    ids = [];
    return;
end

if iscell(namelist{1})
    ids = {};
    for j = 1:length(namelist)
        [names{j}, ids{j}] = FindExpts(namelist{j},varargin{:});
    end
    return;
end

findstr = {};
ids = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'test',4)
    elseif ischar(varargin{j})
        findstr = {findstr{:} varargin{j}};
    end
    j = j+1;
end
nf = 0;
for k = 1:length(findstr)
for j = 1:length(namelist)
    if strfind(namelist{j},findstr{k})
        nf = nf+1;
        names{nf} = namelist{j};
        ids(nf) = j;
    end
end
end



function [names, details] = ListExptCells(E, varargin)

for j = 1:length(E)
    names{j} = E{j}.Header.expname;
    if isfield(E{j},'Trials')
        details.ntrials(j) = length(E{j}.Trials);
    elseif isfield(E{j},'ntrials')
        details.ntrials(j) = E{j}.ntrials;
    end
    fprintf('%d: (%dTrials) %s\n',j,details.ntrials(j),names{j});
end

function [names, details] = ListExptDir(E, varargin)

d = dir([E '/*idx.mat']);
nx = 0;
names = {};
details = [];
for j = 1:length(d)
    name = [E '/' d(j).name];
    clear ExptList;
    load(name);
    if ~exist('ExptList','var')
        ExptList = BuildExptList(Expt, Expts);
    end
    for k = 1:length(ExptList)
        e = find([Expts.start] == ExptList(k).start);
        if isempty(Expts(e).result)
            Expts(e).result = 2;
        end
        if length(e) && Expts(e).result == 2
        nx = nx+1;
        names{nx} = ExptList(k).expname;
        details.ntrials(nx) = Expts(e).lasttrial-Expts(e).firsttrial;
        details.filename{nx} = d(j).name;
        end
    end
end


function ExptList = BuildExptList(Expt, Expts)

ExptList = [];

nx = 0;
for j = 1:length(Expts)
    if Expts(j).lasttrial - Expts(j).firsttrial > 10
        nx = nx+1;
        ts = Expts(j).firsttrial;
        ExptList(nx).et = Expt.Trials.et{ts};
        ExptList(nx).e2 = Expt.Trials.e2{ts};
        ExptList(nx).e3 = Expt.Trials.e3{ts};
        ExptList(nx).start = Expts(j).start;
        ExptList(nx).expname = [ExptList(nx).et 'X' ExptList(nx).e2 'X' ExptList(nx).e3];
        if strfind(Expt.Trials.OptionCode{ts},'+fS')
            ExptList(nx).expname = [ExptList(nx).expname 'FS'];
        end
    end
end

function [names, details] = ListExptDirs(E, varargin)

names = {};
details = {};
d = dir(E);
if isdir(E)
    root = [E '/'];
else
    id = strfind(E,'/');
    if length(id)
        root = E(1:id(end));
    end
end
nd = 0;
for j = 1:length(d)
    if d(j).isdir
        nd = nd+1;
        path = [root d(j).name];
        [names{nd}, details{nd}] = ListExptDir(path, varargin);
        details{nd}.dirpath = path;
    end
end
