function [name, details] = GetName(C, varargin)
% Getname(X) returns the cell name associated with string/structre X
% e.g. GetName(Expt) will return lemM001 if thats the correct Expt session
%Getname(X,'path') returns a full pathname
%Getname(X,'folder') returns a the folder the file is in
%if X is a cell array, builds a cell array of names.  If
% this has one unique value, returns this value. Otherwise returns the cell
% array
%GetName(Expt,'withcell')  includes cell number
%GetName(Expt,'withexpt') includes experiment type label
%GetName(Expt,'withsuffix') includes suffix # 
%
% see also BuildFileName
fullpath = 0;
cellid = 0;
showexpt = 0;
nametype = '';
silent = '-default';
loadname = '';
showsuffix = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'path',4)
        fullpath = 1;
    elseif strncmpi(varargin{j},'expt',4)
        nametype = 'expt';
    elseif strncmpi(varargin{j},'folder',6)
        fullpath = 2;
    elseif strncmpi(varargin{j},'lfp',3)
        nametype = 'lfpfile';
    elseif strncmpi(varargin{j},'-silent',3)
        silent = varargin{j};
    elseif strncmpi(varargin{j},'withexpt',7)
        showexpt = 1;
    elseif strncmpi(varargin{j},'withsuffix',7)
        showsuffix = 1;
    elseif strncmpi(varargin{j},'withcell',7)
        nametype = 'withcell';
    elseif strncmpi(varargin{j},'loadname',7)
        nametype = 'loadname';
    end
    j = j+1;
end


cellid = GetCellNumber(C);
details = [];
name = '';
exptno = 0;
if isstruct(C)
    if isfield(C,'Expt') %Expt
        [name, details] = GetName(C.Expt, varargin{:});
    elseif isfield(C,'Header') %Expt
        ename = GetEval(C,'name');
        if fullpath || strcmp(nametype,'expt')
            name = ename;
        elseif strcmp(nametype,'loadname')
            if isfield(C.Header,'loadname')
                name = C.Header.loadname;
            else
                name = ['??' ename];
            end
        else
            [a,b,c,d] = GetMonkeyName(ename);
            name = [a c];
        end
        exptno = GetExptNumber(C);
    elseif isfield(C,'dirname') %RF fit
        [a,b,c,d] = GetMonkeyName(C.dirname);
        name = [a c];
    elseif isfield(C,'V') %FullV file
        name = '';
        if isfield(C,'loadname')
            [a,b,c,d] = GetMonkeyName(C.loadname);
            loadname = C.loadname;           
        elseif isfield(C,'prefix')
            [a,b,c,d] = GetMonkeyName(C.prefix);
            name = [C.prefix '/'];
        elseif isfield(C,'name')
            [a,b,c,d] = GetMonkeyName(C.name);
        end
        if fullpath
            if isempty(name)
                name = [GetFilePath('data') '/' d];
            end                
        else
            name = [a c];
        end
    elseif isfield(C,'loadname') %could be many thing. including RF fix list
        if strcmp(nametype,'loadname')
            name = C.loadname;
        else
            [a,b,c,d] = GetMonkeyName(C.loadname);
            name = [a c];
        end
        exptno = GetExptNumber(C.loadname);
        loadname = C.loadname;
    elseif isfield(C,'name') %could be many thing. including RF fix list
        loadname = C.name;
        [a,b,c,d] = GetMonkeyName(C.name);
        name = [a c];
    elseif isfield(C, 'spkfile') % a cluster without loadname
        [a,b,c,d] = GetMonkeyName(C.spkfile);
        name = [a c];
        bname = regexprep(C.spkfile,'.*/data/','');
        bname = regexprep(bname,'/Spikes/.*','');
        exptno = GetExptNumber(C.spkfile);
        loadname = sprintf('/b/data/%s/Expt%dClusterTimes.mat',bname,exptno);
    end
elseif ischar(C)
        [a,b,c,d] = GetMonkeyName(C, silent);
        name = [a c];
elseif iscell(C)
    for j = 1:length(C)
        names{j} = GetName(C{j},varargin{:});
    end
    uname = unique(names);
    if isempty(uname{1})
        uname = uname(2:end);
    end
    if length(uname) ==1
        name = uname{1};
        details.names = names;
    else
        name = names;
    end
    return;
end
if showsuffix
    name = sprintf('%s.%d',name,exptno);
end

if strcmp(nametype,'loadname') && ~isempty(loadname)
    name = loadname;
end

if fullpath ==2
    if ~isempty(loadname)
        name = fileparts(loadname);
    else
        name = fileparts(name);
    end
end


if strcmp(nametype,'lfpfile')
    eid = GetExptNumber(C);
    name = [name '.' num2str(eid) '.lfp.mat'];
elseif strcmp(nametype,'withcell')
    name = [name '.cell' num2str(cellid)];
end

if showexpt 
    if ischar(C)
        estr = regexprep(d, ['.*' name '.'],'');
    else
        estr = Expt2Name(C);
    end
        name = [name '.' estr];
end

