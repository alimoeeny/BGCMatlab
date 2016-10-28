function str = BuildFileName(name, type, varargin) 
%str = BuildFileName(name, type, ...) build name for an object.
% returns string naming file given type
% types are 'fullv','combine','rffits' 'datadir' 'celllist' 
% 'error'  File containing Seesion Errors
% 'spkfile'
% 'rffile'
%'expt'
%'allexpt'
%'shortname'  Short string identifying session
%'paramfile' File with manually specified parameters (area, electrode)
%'spike2mat'  raw .mat file made by spike2
%'spike2idx'  idx.mat file made by spike2
%'rawdir'   Folder with .smrc files (may be different from .mat when online
str = [];
Expt = [];
AllExpt = [];
probe  =0;
sudir = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'probe',5)
        j = j+1;
        probe = varargin{j};
    elseif strncmpi(varargin{j},'su',2)
        sudir = 1;
    end
    j = j+1;
end
if iscell(name)
    for j = 1:length(name)
        str{j} = BuildFileName(name{j},type, varargin{:});
    end
    return;
end


eid = GetExptNumber(name);
if ~ischar(name)
    S = name;
    dataroot = GetName(S, 'loadname');
    if ~isdir(dataroot)
        dataroot = fileparts(dataroot);
    end
    prefix = GetName(S);
    eid = GetExptNumber(S);
else
    dataroot = fileparts(name);
end
if isfield(name,'Expt') && isfield(name.Expt,'Header')
    AllExpt = name;
    name = GetEval(name.Expt,'name');
elseif isfield(name,'Header')
    Expt = name;
    name = GetEval(name,'name');
elseif ~ischar(name)
    if strcmp(type,'error')
        name = [GetName(name,'folder') '/'];
    else
        name = GetName(name);
    end
end

name = regexprep(name,'\\','/');
if nargin == 1
    type = 'datadir';
end
   
drive = [];
if regexp(name,'^[A-Z]:')
    drive = name(1:2);
    if isunix
        if strcmp(name(4:7),'data')
            newname = ['/b/data' name(8:end)];            
            fprintf('BuildFileName: Replacing %s/data with /b/data for Mac/Unix',drive);
            name = newname;
        end
    end
end

[a,b,c,d] = GetMonkeyName(name);
if isdir(name)
    name = [name '/']; 
end
[root, file] = fileparts(name);



if strcmp(type,'fullv')
    if probe
        str = [fileparts(name) '/Expt' num2str(eid) '.p' num2str(probe) 'FullV.mat'];
    else
        str = [fileparts(name) '/Expt' num2str(eid) 'FullV.mat'];
    end
elseif strcmp(type,'shortname')
    str = [b c];
elseif strcmp(type,'paramfile')
    str = [root '/' b c '.params'];
elseif strcmp(type,'autocluster') % make filename to give to combine
    [root, file] = fileparts(name);
    if isnumeric(varargin{1})
        str = [root '/Expt' int2str(varargin{1}) 'AutoClusterTimes.mat'];
    end
elseif strcmp(type,'expt')
    if ~isempty(Expt)
        str = [root '/' a c '.mat'];
    end
elseif strcmp(type,'allexpt')
    if ~isempty(AllExpt)
        str = [root '/' a c '.' Expt2Name(AllExpt.Expt) '.Cells.mat'];
    end
elseif strcmp(type,'combine') % make filename to give to combine
    [root, file] = fileparts(name);
    if c(1) == 'M' %laminar probe - use directory
        str = root;
    else
        str = [root '/' a c '.mat'];
    end
elseif strcmp(type,'datadir')
    str = regexprep(name,['/' c '/.*'],['/' c]);
    if isempty(strfind(name,'/'))
        if sudir
            str = [drive '/b/data/' a '/SU/' c '/'];
        else
            str = [drive '/b/data/' a '/' c '/'];
        end
    end
elseif strncmp(type,'error',5)
    str = [root '/Errors.mat'];
elseif strcmp(type,'celllist')
    str = regexprep(name,['/' c '/.*'],['/' c]);
    if isempty(strfind(name,'/'))
        str = ['/b/data/' a '/' c '/'];
    end
    str = [str '/CellList.mat'];
elseif strcmp(type,'rffits')
    if regexp(name,[a '/SE/'])
        dataroot = regexprep(name,['data/' a '/SE/.*' ],['data/' a '/SE/']);
    else %was once c instead of name. Surely wrong = e.g. /b/data/jbe/M* c is M*
        dataroot = regexprep(name,['data/' a '/.*' ],['data/' a]);
    end
    str = [dataroot '/rffits.mat'];
elseif strcmp(type,'bnc')
    str = regexprep(name,['/' c '/.*'],['/' c]);
elseif strcmp(type,'rffile')
    if regexp(name,[a '/SE/'])
        dataroot = regexprep(name,['data/' a '/SE/.*' ],['data/' a '/SE/']);
    else
        dataroot = regexprep(name,['data/' a '/.*' ],['data/' a]);
    end
    str = [dataroot '/' c '/' a c '.rf.mat'];
elseif strcmp(type,'spkfile')
    dataroot = GetName(S, 'loadname');
    prefix = GetName(S);
    if ~isdir(dataroot)
        dataroot = fileparts(dataroot);
    end
    p = GetProbeNumber(S);
    e = GetExptNumber(S);
    str = sprintf('%s/Spikes/%s.p%dt%d.mat',dataroot,prefix,p,e);
    if ~exist(str) %?try differen path. Not if using this to build name. ?use flag if want to check
    end
elseif strcmpi(type,'spikes')
    str = sprintf('%s/Expt%dSpikes',dataroot,eid);    
elseif strcmpi(type,'exptsummary')
    str = sprintf('%s/ExptSummary.mat',root);
elseif strcmp(type,'rawdir')
    str = sprintf('%s/Spike2/data/%s/%s',dataroot(1:2),a,c);
elseif strcmp(type,'spike2mat')
    str = sprintf('%s/%s.%d.mat',dataroot,prefix,eid);
elseif strcmp(type,'spike2idx')
    str = sprintf('%s/%s.%didx.mat',dataroot,prefix,eid);
end
