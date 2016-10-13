function details = listop(list, op, varargin)
%
% listop(list, 'copy', ....)
% copies all files in list, on in the current drive to C:
% listop(list, 'copy', 'cpem',....)
%   copies over any .em.mat files that match files in list
%
%
% copies files
% Do useful things to files named in list.

cprefix = 'C:';
skip = 1;
doac = 1;
doorfiles = 1;
doMTfiles = 0;
doemfiles = 0;
funcfcn = [];
funcargs = {};
tgtdir = 'C:/data';

if iscellstr(list)
    fstrings = list;
elseif ischar(list)
    if ~exist(list)
        fprintf('Cant Read %s\n',list);
        return;
    else
        fstrings = textread(list,'%s','delimiter','\n');
    end
end
Expts = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'sizestrings',7)
        OTTF.sizestrings = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'cpem',4)
        doemfiles = 1;
    elseif strncmpi(varargin{j},'function',3)
        j = j+1;
        funcfcn = varargin{j};
        if length(varargin) > j & iscell(varargin{j+1})
            j = j+1;
            funcargs = varargin{j};
        end
    elseif strncmpi(varargin{j},'cylinder',5)
        doMTfiles = 1;
        doorfiles = 0;
        doacfiles = 0;
    elseif strncmpi(varargin{j},'size',2)
        j = j+1;
        wsize = varargin{j};
    elseif strncmpi(varargin{j},'skip',2)
        j = j+1;
        skip = varargin{j};
    elseif strncmpi(varargin{j},'tgtdir',6)
        j = j+1;
        tgtdir = varargin{j};
    end
    j = j+1;
end

details.ndone = 0;
ndone = 0;
names = [];

if strncmpi(op,'removemu',8)
        details = removemufiles(fstrings);
        return;
end

rid = strmatch('#',fstrings);
id = setdiff(1:length(fstrings),rid);
fstrings = fstrings(id);


for j = skip:length(fstrings)
if length(fstrings{j}) == 0;
elseif length(funcfcn)
    pathname = fstrings{j};
    if isempty(funcargs)
        good = eval([funcfcn '(''' pathname ''')']);
    else
        str = [];
        for k = 1:length(funcargs)
            str = [str ',''' funcargs{k} ''''];
        end
        good = eval([funcfcn '(''' pathname '''' str ')']);
    end
elseif strncmpi(op,'recombineorbw',12)
    name = deblank(regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat'));
    name = regexprep(name,'\.cell[0-9]*\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) % not a brainwave file
        combine(name,'recombineorbw','recount');
    end
elseif strncmpi(op,'allprobes',6)
    if regexp(fstrings{j},'\.cell[0-9]*\..*.mat')
    name = regexprep(fstrings{j},'\.cell[0-9]*\..*.mat','.mat');
    if ndone == 0  | isempty(strmatch(name,names))
        PlotAllProbes(fileparts(name),'save');
    ndone = ndone+1;
    names{ndone} = name;
    end
    end
    details.names = names;
    details.ndone = ndone;
elseif strncmpi(op,'load',4)
    ndone = ndone+1;
    Expts{ndone} = LoadExpt(fstrings{j});
    details.Names{ndone} = fstrings{j};    
    if expt.isall(Expts{ndone});
        A = All2Expt(Expts{ndone});
        for k = 1:length(A)
            Expts{ndone} = A{k};
            Expt = a{k};
            details.Names{ndone} = GetName(A{k},'withcell');
            x = whos('Expt');
            details.size(ndone) = x.bytes;
            ndone = ndone+1;
        end
        ndone = ndone-1;
    else
        Expt = Expts{ndone};
         x = whos('Expt');
         details.size(ndone) = x.bytes;
    end
elseif strncmpi(op,'Check',4)
    ndone = ndone+1;
    fprintf('Reading %d %s',j,fstrings{j}); 
    Expt = LoadExpt(fstrings{j});
    x = whos('Expt');
    fprintf(' (%d bytes)\n',x.bytes);
    Expts(ndone).size = x.bytes;
    if isfield(Expt,'Trials')
        Expts(ndone).nt = length(Expt.Trials);
    else
        Expts(ndone).nt = 0;
    end
    if ndone > 1000
        clear Expt;
    end
elseif strncmpi(op,'recombineall',11)
    name = regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) % not a brainwave file
        combine(name,'recombinequit','recount');
    end
elseif strncmpi(op,'recombine',6)
    name = regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) % not a brainwave file
        combine(name,'recombinename',fstrings{j});
    ndone = ndone+1;
    end
elseif strncmpi(op,'chkcluster',6)
    f = load(fstrings{j});
    if isfield(f,'cExpt')
        f.cExpt.Header.Filename = fstrings{j};
        [dp(j), d] = ExptCellQuality(f.cExpt,'verbose');
        ndone = ndone+1;
    elseif isfield(f,'Expt')
        Expt.Header.Filename = fstrings{j};
        [dp(j), d] = ExptCellQuality(Expt,'verbose');
        ndone = ndone+1;
    end
    details.errs(j) = d.errs;
        
elseif strncmpi(op,'relist',6)
    name = regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat');
    name = regexprep(name,'\.cell[0-9]\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) & (ndone == 0  | isempty(strmatch(name,names)))
        APlaySpkFile(name,'relist');
        ndone = ndone+1;
        names{ndone} = name;
    end
    details.names = names;
    details.ndone = ndone;
elseif strncmpi(op,'copy',4)
    tgt = [cprefix fstrings{j}];

    a = CopyFile(fstrings{j},tgt);
    ndone = ndone+a;
    details.ndone = ndone;
    if doemfiles
        emfile = regexprep(fstrings{j},'\..*.mat','.em.mat');
        if ~strcmp(emfile,fstrings{j}) & exist(emfile,'file');
        tgt = [cprefix emfile];
        CopyFile(emfile,tgt);
        end
    end
    if doorfiles
        otfile = regexprep(fstrings{j},'image.ORBW.','image.OT.');
        if ~strcmp(otfile,fstrings{j}) & exist(otfile,'file');
        tgt = [cprefix otfile];
        CopyFile(otfile,tgt);
        end
        otfile = regexprep(fstrings{j},'image.ORBW.','grating.OXM.');
        if ~strcmp(otfile,fstrings{j}) & exist(otfile,'file');
        tgt = [cprefix otfile];
        CopyFile(otfile,tgt);
        end
    end
    if doMTfiles
        id = regexp(fstrings{j},'[0-9,A-Z][A-Z][A-Z]*\.mat');
        suffix = fstrings{j}(id:end-4);
        acfile = regexprep(fstrings{j},suffix,'DT');
        if ~strcmp(acfile,fstrings{j})
            CopyFile(acfile,[cprefix acfile]);
        end
        acfile = regexprep(fstrings{j},'cylinder.DID','rds.OT');
        if ~strcmp(acfile,fstrings{j})
            CopyFile(acfile,[cprefix acfile]);
        end
    end
    if doac
        id = regexp(fstrings{j},'[0-9,F][A-Z][A-Z]*\.mat');
        suffix = fstrings{j}(id:end-4);
        if isempty(id)
            id = regexp(fstrings{j},'[0-9,F][A-Z][A-Z]*\.cell[0-9]*.mat');
            if isempty(id)
                id = regexp(fstrings{j},'[0-9,F][A-Z][A-Z]*\.mu[0-9]*.mat');
                cid = regexp(fstrings{j},'\.mu[0-9]*.mat');
            else
                cid = regexp(fstrings{j},'\.cell[0-9]*.mat');
            end
            if ~isempty(id)
                suffix = fstrings{j}(id:cid(1)-1);
            else
                suffix = [];
            end
        end
        if ~isempty(suffix)
        [a,b] = regexp(fstrings{j},'\.[a-z][a-z]*\.','start','end');
        stim = fstrings{j}(a+1:b-1);
        acfile = regexprep(fstrings{j},suffix,'AC');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},suffix,'OXAC');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},['\.' stim '\.' suffix],'.nsines.DP');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},['\.' stim '\.' suffix],'.rds.AC');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},['\.' stim '\.' suffix],'.rds.DT');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        end
    end
elseif doemfiles
   emfile = regexprep(fstrings{j},'\.c[0-9]\.','.em.');
   if emfile(1) ~= '/'
       [mnk, monkey, pref, b] = GetMonkeyName(fstrings{j});
       a = strrep(b,mnk,[monkey '/']);
       src = ['Y:/data/' a '/' emfile ];
   end
   tgt = [tgtdir '/' emfile];
   [a,b,c] = copyfile(src, tgt);
   if a > 0 
       ndone = ndone+1;
   end
end
end
if length(Expts)
    details.Expts = Expts;
end
details.ndone = ndone;


function newlist = removemufiles(list)


su = ones(size(list));
for j = 1:length(list)
    if regexp(list{j},'\.p[0-9]+c[0-9]+')
        su(j) = 0;
    elseif regexp(list{j},'.lfp.')
        su(j) = 0;
    elseif regexp(list{j},'\.em\.')
        su(j) = 0;
    elseif regexp(list{j},'\.mu[0-9]+\.')
        su(j) = 0;
    elseif regexp(list{j},'backup')
        su(j) = 0;
    end
    
end
newlist = list(find(su));

function done = CopyFile(src, tgt, varargin)
%copy src to tgt, if tgt doesn't exist or is older than src.
done = 0;
sd = dir(src);
if strcmp(src,tgt)
    fprintf('%s and %s are the same\n',src, tgt);
    return;
end    
if ~exist(src)
    fprintf('Cant Find %s\n',src);
    return;
end
if isempty(sd)
    fprintf('Cant Find %s\n',src);
    return;
end
    d = dir(tgt);
    dd = dir(fileparts(tgt));
    if length(d) == 0  || d(1).datenum < sd(1).datenum
        if length(dd) == 0
           mkdir(fileparts(tgt));
        end
        copyfile(src, tgt);
        cmd = ['copy ' src ' ' tgt];
        cmd = strrep(cmd,'/','\');
        fprintf('%s\n',cmd);
        done = 1;
    end

