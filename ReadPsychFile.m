function [Expts, Summary] = ReadPsychFile(name, varargin)
%[Expts, Summary] = ReadPsychFile(name, varargin)
%   .   ,'useallexpts');  makes expts for psych and non-psych blocks
%   .... ,'nmin',n) sets min # trials to make expt block
%
%ReadPsychFile(name, 'noexpts') returns a signle data structure with all
%Trials
%
%File format
% R%d  0 = wrong
%      1 = correct
%      2 = late/foul
%      3 = BadFix
%      >4 special lines with stimulus/expt info
%      5 start block
%      4 ednd block
%      9 start block of Human psych
%      11 Badfix caused by a microsaccade
%      51 completed fixation trial
%      53 badfix Trial

% R%d xx=%f yy=%f time trialdur rwsize
% xx = value for expt 1
% yy = value for expt 2
%
% R = 0 = WRONG, 1 = CORRECT, 2 = FOUL Choice 3 = BAD-FIX
% R = 100 or 101 are 0,1 but in a correction loop
% 50,51 means saccade was not required
% 53 = fixation only trial, but bad fixation
% R=109 indicates start of expt
% R9 Expt Start for Human Psych 
% R4 stimulus properties at start of expt, 
% R5 other stimulus properties
% R7 stimulus properties not in strict format - don't send these lines to
%                         textscan
% R27 Cancel Exp
% R10 = Expt End
% R8 = Expt Finished by verg


buildexpts = 1;
startday = 0;
DATA.verbose = 0;
DATA.plot.round = [0 0]; 
DATA.filename = name;
DATA.plot.mintrials = 10;
DATA.useallexpts = 0;
fixfield = '';
Expts = {};
Summary.ntrials = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noexpts',8)
        buildexpts = 0;
    elseif strncmpi(varargin{j},'fix',3)
        j = j+1;
        fixfield = varargin{j};
    elseif strncmpi(varargin{j},'useallexpts',8)
        DATA.useallexpts = 1;
    elseif strncmpi(varargin{j},'nmin',4)
        j = j+1;
        DATA.plot.mintrials = varargin{j};
    end
    j = j+1;
end

if iscellstr(name)
    for j = 1:length(name)
        fprintf('%d:reading %s\n',j,name{j});
        Expts{j} = ReadPsychFile(name{j},varargin{:});
    end
    return;
end
tic;
txt = scanlines(name);
d = dir(name);
if isempty(d)
    return;
end
CreationDate = d.datenum; %in case can't find date in file
gid = find(strncmp('R0',txt,2) | strncmp('R1',txt,2) | strncmp('R4',txt,2) | strncmp('R5',txt,2));
xid = setdiff(1:length(txt),gid);
gid = find(strncmp('R',txt,1));

if isempty(gid)
    cprintf('red','No Data in %s\n',name);
    return;
end

DATA.gid = gid;

if ~isempty(xid)
    vstr = 've=';
    ids = strfind(txt(xid),vstr);
    gotve = find(CellToMat(ids));
    
    if isempty(gotve)
        vstr = 'binoclean=';
        ids = strfind(txt(xid),vstr);
        gotve = find(CellToMat(ids));
    end
    if isempty(gotve)
        DATA.binocversion(1) = 0;
    else
    ve = sscanf(txt{xid(gotve(1))}(ids{gotve(1)}(1):end),[vstr '%f.%f']);
    DATA.binocversion(1) = ve(1);
    if length(ve) > 1
        DATA.binocversion(2)= ve(2);
    end
    end
    vstr = 'date=';
    dids = find(CellToMat(strfind(txt(xid),vstr)));
    if ~isempty(dids)
    id = findstr(txt{xid(dids(1))},'date=');
    DATA.recorddate = datenum(txt{xid(dids(1))}(id+9:id+28));
    end
else
    DATA.binocversion(1) = 0;
end
if ~isfield(DATA,'recorddate')
    [~,s] = fileparts(name);
    DATA.recorddate = datenum(s(4:end));
end
Stimulus.rw = 0;
Stimulus.RespDir = 0;
Stimulus.tid = 0;
Stimulus.rwdir = 0;
nt = 0;
nx = 0;
Expts(1).start = 0;
ptime = 0;
bwtime = 0;
txt = regexprep(txt,'(^R[0-9]+)S ','$1 ');
tfield = zeros(1, length(gid));
for j = 1:length(gid)
    atxt = txt{gid(j)};
    txt{gid(j)} = regexprep(txt{gid(j)},'mo=[a-z]+ ','');
    x = textscan(txt{gid(j)}, 'R%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f');
    if ~isempty(x{8}) && ischar(x{8}{1}) && ~isempty(regexp(x{8}{1},'^[0-9]+'))
        ttimes = sscanf(x{8}{1},'%f %f %f');
        x{8} = regexprep(x{8},'.* ','');
        tfield(j) = 8;
    end
%    fprintf('%d%s\n',j,atxt);
    DATA.score(j) = x{1}(1);
    DATA.times(j) = 0;
    s = txt{gid(j)};
    if DATA.score(j) == 7
        id = strfind(s,'date=');
        if ~isempty(id)
            CreationDate = datenum(s(id(1)+[9:23]));
        end
        id = strfind(s,'progtime=');
        if ~isempty(id)
            ptime = sscanf(s(id(1):end),'progtime=%f');
        end
        id = strfind(s,'bt=');
        if ~isempty(id)
            bwtime = sscanf(s(id(1):end),'bt=%f');
        end
    else
        if length(x) < 12
        end
        if DATA.score(j) == 4
        end
    for k = 2:length(x)        
        if ~isempty(x{k})
        if (k > 2 && mod(k,2) == 0 && ~isempty(regexp(x{k}{1},'^[0-9]'))) ...
            || tfield(j) == k
            if tfield(j) == 0
                tfield(j) = k;
                val = sscanf(x{k}{1},'%f');
                goodf = 0;
            else %three floats AND the next field
                val = ttimes;
                goodf = 1;
            end
            if DATA.score(j) == 4
                nx = nx+1;
                Expts(nx).start = val(1);
                Expts(nx).firsttrial = nt+1;
            elseif DATA.score(j) == 5 
                if nx > 0
                Expts(nx).End = val(1);
                Expts(nx).lasttrial = nt;
                end
            else
                DATA.times(j) = val(1);
            end
%if is quite a bit faster....            
            if DATA.score(j) ==1 || DATA.score(j) == 0 || DATA.score(j) == 51 || DATA.score(j) == 53
%           if ismember(DATA.score(j),[0 1 51 53])
                DATA.rwszs(j) = val(3);
                if val(3) < 0.001 || val(3) > 0.8
                    val(3) = 0;
                end
            else
                DATA.rwszs(j) = 0;                
            end
            DATA.DURS(j) = val(2);
        elseif k > 2 && mod(k,2) == 0 && strcmp(x{k}{1},'sn')
            DATA.sign(j) = x{k+1};
            Stimulus.rwdir = DATA.sign(j);
            goodf = 1;
        else
            goodf = 1;
        end
        if iscell(x{k})
            a{k}{j} = x{k}{1};
            if goodf
                lastf = x{k}{1};
            else
                lastf = '';
            end
            if strcmp(lastf,'L')
                Stimulus.RespDir = -1;
                lastf = '';
            elseif strcmp(lastf,'R')
                Stimulus.RespDir = 1;
                lastf = '';
            elseif strncmp(lastf,'Sa',2)
                lastf = 'Sa';
            end
        else
            if isvarname(lastf)
                Stimulus.(lastf) = x{k};
            end
            lastf = '';
            a{k}(j) = x{k};
        end    
        end
    end
    end
    Stimulus.saved = strncmp('S ',x{2},2);
    if ismember(DATA.score(j),[0 1 50 51])
        snid = strfind(txt{gid(j)},'stid=');
        if ~isempty(snid)
            Stimulus.stid = sscanf(txt{gid(j)}(snid+5:end),'%f');
        end
            
        if length(tfield) < j
            tfield(j) = 0;
        end
        if Stimulus.rwdir == 0 
            snid = strfind(txt{gid(j)},'sn=');
            if ~isempty(snid)
                rwd = sscanf(txt{gid(j)}(snid+3:end),'%f');
                Stimulus.rwdir = rwd(1);
            end            
        end
        
        
        nt = nt+1;
        snid = regexp(txt{gid(j)},' [0-9]+[0-9,.]* ');
        if ~isempty(snid)
            vals = sscanf(txt{gid(j)}(snid+1:end),'%f');
            Stimulus.Start = vals(1)-ptime;
            DATA.times(j) = vals(1);
            if length(vals) > 1
                Stimulus.End = vals(1)+vals(2)-ptime;
            end
        end
        if Stimulus.RespDir == 0 && Stimulus.rwdir ~= 0
            resp = (DATA.score(j) .* 2)-1;
            Stimulus.RespDir = Stimulus.rwdir .* resp;
        end
        Trials{nt} = Stimulus;
        Trials{nt}.xtype = x{2}{1};
        Trials{nt}.ytype = x{4}{1};
        Trials{nt}.linenumber = j;
        Trials{nt}.date = CreationDate + (DATA.times(j) - ptime)./(24 * 60 * 60);
        if nt > 1 && Stimulus.tid > 0 && Stimulus.tid == Trials{nt-1}.tid
            Trials{nt-1}.RespDir = Trials{nt-1}.RespDir .*2;
        end
        if tfield(j) > 6
            Trials{nt}.ztype = x{6}{1};
        elseif tfield(j) == 0
            Trials{nt}.ztype = 'e0';
        end
        Stimulus.Sa = 0;
        Stimulus.RespDir = 0;
        Stimulus.rwdir = 0;
        Stimulus.id = 0;
        Stimulus.Start = 0;
        Stimulus.End = 0;
    end
end    
if fixfield
    Trials = FixField(Trials,fixfield);
end
if exist('Trials','var');
    DATA.Trials = CellToStruct(Trials);
    DATA.saved = [DATA.Trials.saved];
end

id = regexp(name,'[0-9][0-9][A-Z]');
if length(id)
    startday = datenum(name(id(1):id(1)+8));
end
if DATA.verbose
    fprintf('\n%d lines,', length(a{1}));
end
if isnumeric(a) & a == -1 %textscan returned an error
    DATA.score = [];
    return;
end
if isempty(a{8})
    return;
end
d = dir(name);
DATA.filedate = d.datenum;
nxval = (length(a)-9)/2;

%
% find experiment actual date
id = find(CellToMat(strfind(txt(xid),'date=')));
progtimeoffset = 0;
ptime = 0;
bwtime = 0;
if ~isempty(id)
    s = txt{xid(id(1))};
    id = strfind(s,'date=');
    DATA.CreationDate = datenum(s(id(1)+[9:23]));
    id = strfind(s,'progtime=');
    if ~isempty(id)
        ptime = sscanf(s(id(1):end),'progtime=%f');
    end
    id = strfind(s,'bt=');
    if ~isempty(id)
        bwtime = sscanf(s(id(1):end),'bt=%f');
    end
    if ptime > 0 && bwtime > 0
        progtimeoffset = ptime-bwtime;
    end
    
else
    id = find(CellToMat(strfind(txt(xid),'time=')));
    if ~isempty(id)
        s = txt{xid(id(1))};
        id = strfind(s,'time=');
        h = sscanf(s(id(1):end),'time=%d:%d');
        DATA.CreationDate = startday + (h(1) + h(2)./60)./24;
    end
end

if length(a) > 10
id = strmatch('Sa',a{11});  %%trials terminated by mirosaccade
if max(id) <= length(DATA.score)
tid = find(DATA.score(id) == 3);
DATA.score(id(tid)) = 7;
end
end
if DATA.verbose
    fprintf('%d xvals,', nxval);
end

[nx, ptypes] = Counts(a{2});
[~,id] = sort(nx);

DATA.xtype = a{2};
if DATA.plot.round(1)
    DATA.x = round(DATA.x/DATA.plot.round(1)) .* DATA.plot.round(1);
end
[nx, ptypes] = Counts(a{4});
[~,id] = sort(nx);
DATA.ytype = a{4};
if DATA.plot.round(2)
    DATA.y = round(DATA.y/DATA.plot.round(2)) .* DATA.plot.round(2);
end



id = find(strcmp('e0',a{4}));
DATA.y(id) = 0;
%DATA.times(find(DATA.score == 4)) = 0;
if length(a) > 11
    DATA.ztype = a{11};
    DATA.z = a{12};
end

if DATA.times(1) == 0
    id = find(DATA.times > 0 & DATA.score ~= 6);
    DATA.times(1) = DATA.times(id(1));
end
esid = find(a{1} ==4);
for j = 1:length(esid)
    if a{8}(esid(j)) > 1000
        DATA.times(esid(j)) = ctime2datenum(a{8}(esid(j)));
    end
end

id = find(DATA.times ==0);
while ~isempty(id)
    DATA.times(id) = DATA.times(id-1);
    id = find(DATA.times ==0);
end
DATA.times = DATA.times-progtimeoffset;
    
if DATA.verbose
    fprintf('%d jumps,', j);
end

DATA.readtime = now;

tid = find(ismember(DATA.score, [4])); %start expts
oid = strmatch('tr',a{4}(tid));
trs = a{5}(tid(oid));
DATA.Blockid.tr = tid(oid);
DATA.Block.tr = trs;
DATA.Stimvals.tr = median(trs);  

  
tid = find(DATA.score == 5); %stim values
oid = strmatch('or',a{2}(tid));
ors = a{3}(tid(oid));
DATA.Blockid.or = tid(oid);
DATA.Block.or = ors;
DATA.Stimvals.or = median(ors);  
oid = strmatch('sz',a{2}(tid));
ors = a{3}(tid(oid));
DATA.Blockid.sz = tid(oid);
DATA.Block.sz = ors;
DATA.Stimvals.sz = mean(ors);  

oid = strmatch('co',a{4}(tid));
cos = a{5}(tid(oid));
DATA.Blockid.co = tid(oid);
DATA.Block.co = cos;
DATA.Stimvals.co = mean(cos);  
oid = strmatch('sf',a{4}(tid));
ors = a{5}(tid(oid));
DATA.Blockid.sf = tid(oid);
DATA.Stimvals.sf = median(ors);  
DATA.Block.sf = ors;

if length(a) > 15
oid = strmatch('dd',a{15}(tid));
if length(oid)
    xs = a{16}(tid(oid));
    DATA.Blockid.dd = tid(oid);
    DATA.Block.dd = xs;
    DATA.Stimvals.dd = mean(xs);  
end
end

if length(a) > 17
oid = strmatch('c2',a{17}(tid));
if length(oid)
    xs = a{18}(tid(oid));
    DATA.Blockid.c2 = tid(oid);
    DATA.Block.c2 = xs;
    DATA.Stimvals.c2 = mean(xs);  
end
end


DATA.Stimvals.wi = round(DATA.Stimvals.sz * 10)/10;
DATA.Stimvals.hi = round(DATA.Stimvals.sz * 10)/10;

oid = find(strncmp('ve',a{2}(tid),2));
stid = tid(oid);

if ~isempty(oid) && DATA.binocversion(1) == 0
    DATA.binocversion = mean(a{3}(stid));
end

if DATA.verbose
    fprintf('%d ves,', length(oid));
end

Summary.ntrials = sum(ismember(DATA.score,[51 1 0]));
rid = find(ismember(DATA.score,[51 1]));
Summary.goodtrials = length(rid);
if isfield(DATA,'rwszs')
    Summary.TotalReward = sum(DATA.rwszs(rid));
end

[sid,eid] = FindBlocks(DATA);

Comments = [];
if length(oid)
        
    opid = find(CellToMat(strfind(txt(xid),'op=')));
    ves = a{3}(tid(oid));
    DATA.Blockid.ve = tid(oid);
    DATA.Block.ve = ves;
    DATA.Stimvals.ve = mean(ves);
    
    OptionCode = regexprep(txt(xid(opid)),'.*op=','op=');
    OptionCode = regexprep(OptionCode,' .*','');
    opid = find(CellToMat(strfind(txt(xid),'op=')));
    id = MatchInd(opid,oid);
    for j = 1:length(id)
        if id(j) > 0
            DATA.Block.OptionCode{id(j)} = OptionCode{j};
            DATA.Blockid.OptionCode(id(j)) = opid(j);
        end
    end
    ns = 0;
    for j = 1:length(sid)
        lines = sid(j):eid(j);
        id = find(xid > gid(sid(j)) & xid < gid(eid(j))); %R7 lines that apply to this expt
        ns(j) = length(id);
        if ns(j) > 400
            id
        end
        for k = 1:length(id)
            s = split(txt{xid(id(k))});
            if sum(strncmp('st=',s,3));
                l = find(strncmp('st=',s,3));
                DATA.Block.st(j) = StimulusName(s{l(1)}(4:end));
                DATA.Blockid.st(j) = gid(sid(j));
            end
        end            
    end
end

id = find(CellToMat(strfind(txt(xid),'cm=')));
for j = id(:)'
    tid = find(gid < xid(j));
    if isempty(tid)
        tid = 1;
    else
        tid = tid(end);
    end
    p = strfind(txt{xid(j)},'cm=');
    Comments(end+1).text = txt{xid(j)}(p:end);
    Comments(end).time = DATA.times(tid);
end


tid = find(DATA.score == 8); %background values
if length(tid)
oid = strmatch('xo',a{2}(tid));
    
end
oid = strmatch('xo',a{2}(tid));
if length(oid)
    xs = a{3}(tid(oid));
    DATA.Blockid.backxo = tid(oid);
    DATA.Block.backxo = xs;
    DATA.Stimvals.backxo = mean(xs);  
end
oid = strmatch('yo',a{4}(tid));
if length(oid)
    xs = a{5}(tid(oid));
    DATA.Blockid.backyo = tid(oid);
    DATA.Block.backyo = xs;
    DATA.Stimvals.backyo = mean(xs);  
end




tid = find(DATA.score == 6); %seed offset for images values
if length(tid)
    DATA.Stimvals.seedoffset = a{12}(tid(end)); %kludge assumes 4x2 stims. Need to handle more carefully one day
end
if DATA.verbose
    fprintf('%d rws,', length(DATA.rwszs));
end
tid = find(DATA.score == 4); %seed offset for images values
if ~isempty(tid)
    ts = 719528.8+[Expts.start]./(60 * 60 * 24);  %convert sec to days
    for j = 1:length(tid)
        id = find(stid < tid(j));
        if length(id)
            DATA.Blockid.Start(id(end)) = tid(j);
            DATA.Block.Start(id(end)) = ts(j);
        end
    end
end
if isfield(DATA.Block,'Start') && max(DATA.Block.Start) > datenum('01-Jan-2000')
    for j = 2:length(DATA.Block.Start)
        if DATA.Block.Start(j) == 0
            DATA.Block.Start(j) = DATA.Block.Start(j-1);
        end
    end
elseif startday > 0
    DATA.Block.Start = startday; 
    DATA.Blockid.Start = 1; 
end

if length(DATA.rwszs) ~=  length(DATA.score)
  n = min([length(DATA.rwszs)  length(DATA.score)]);
  DATA.rwszs = DATA.rwszs(1:n);
  DATA.score = DATA.score(1:n);
end

DATA.rwsum = cumsum(DATA.rwszs .* (ismember(DATA.score,[1 51])));
if DATA.verbose
    fprintf('%d rws,', length(DATA.rwszs));
end

a = diff(find(DATA.score > 100));  %Corr loops
nconsec = 1;
maxconsec = 1;
for j = 1:length(a)
    if a(j) == 1
        nconsec = nconsec+1;
        if nconsec > maxconsec
            maxconsec = nconsec;
        end
    else
        nconsec = 1;
    end
end
if maxconsec > 20
    fprintf('Max corr loop length %d\n',maxconsec);
end

ids = find(ismember(DATA.score,[0 1 2 3 7 100 101 102 103]));
if length(ids) < length(DATA.score)/5  %not psychophysics
    DATA.ispsych = 0;
else
    DATA.ispsych = 1;
end
if DATA.binocversion == 0
    DATA.Stimvals.ve = 0;
end
if buildexpts
    DATA.Comments = Comments;
    Expts = MakeExpts(DATA,Expts);
else
    Expts = DATA;
end

function Trials = FixField(Trials, field, varargin)
    stid = [];
    for j = 1:length(Trials)
        if isfield(Trials{j},field)
            good(j) = 1;
        else
            good(j) = 0;
        end
        if isfield(Trials{j},'stid')
            stids(j) = Trials{j}.stid;
        end
    end
    if strcmp(field,'Adapt') && ~isempty(stids)
        splitid = max(stids)/2;
        a = ((stids>splitid)-0.5)/5; %+- 0.1
        for j = 1:length(Trials)
            if good(j) ==0
                Trials{j}.Adapt = a(j);
            end
        end
        
    end
    
    
function [sid, eid] = FindBlocks(DATA)
    
    sid = [];
    eid = [];
    if ~isfield(DATA,'score')
    return;
    end

if DATA.binocversion(1) < 0.2 || DATA.binocversion(1) > 4
    id = find(ismember(DATA.score,[4 109]));  %block starts.
    eid = find(DATA.score(1:end-1) == 8 & DATA.score(2:end) == 5); %Expt End or Cancel
else
    id = find(ismember(DATA.score,[9 109]));  %block starts.
    %Get R8 when Stim2Psych is called with flag 0. When Verg forces end
    
    eid = find(ismember(DATA.score,[8 10 27])); %Expt End or Cancel
end
if isempty(eid)
    eid = length(DATA.score);
end



if isempty(id)
    fprintf('No Expts found\n');
    return;
end

while eid(1) < id(1)
    eid = eid(2:end);
end
blkstart = id;
blkend = eid;

if isempty(id)
    id(1) = 1;
    eid(1) = length(DATA.score)-1;
elseif eid(end) < id(end)
    eid(end+1:length(id)) = length(DATA.score)-1;
end
nx = 1;
if length(eid) > length(id) && eid(1) < id(1)
    eid = eid(2:end);
end

j = 2;
while j < length(id)
    if length(eid) >= j-1
        if id(j) < eid(j-1)
            fprintf('Blocks misaligned at line %d\n',DATA.gid(id(j)));
            ng = sum(ismember(DATA.score(id(j-1):id(j)),[0 1 51]));
            if ng < 3
                fprintf('Ignoring Empty Expt Lines  %d - %d\n',DATA.gid(id(j-1)),DATA.gid(id(j)));
                id(j-1) = [];
                j = j-1;
            else
                fprintf('Cant ingore lines %d-%d - has %d trials\n',DATA.gid(id(j-1)),DATA.gid(id(j)),ng);
            end
        end
    else
        fprintf('Blocks misaligned (missing end) at %d\n',DATA.gid(id(j)));        
    end
    j = j+1;
end

sid = id;
if length(eid) < length(sid) %crashes left expts empty
    fprintf('Last end Block at %d, last start at %d\n',DATA.gid(eid(end)),DATA.gid(sid(end)));
    eid(1:length(sid)-1) = sid(2:end);
    eid(length(sid)) = length(DATA.score);
end


function E = FillStimvals(E)
    f = fields(E.Trials);
    for j = 1:length(f)
        vals = [E.Trials.(f{j})];
        if length(unique(vals)) == 1
            E.Stimvals.(f{j}) = unique(vals);
            E.Trials = rmfield(E.Trials,f{j});
        end
    end
        
function Expt = MakeOneExpt(DATA, E)
    if isfield(E,'lasttrial')
        id = E.firsttrial:E.lasttrial;
    else
        id = E.firsttrial:length(DATA.Trials);
    end
    if isempty(id)
        Expt = [];
        return;
    end
    Expt.Trials = DATA.Trials(id);
    Expt.Stimvals.et = cellmax({DATA.Trials(id).xtype});
    Expt.Stimvals.e2 = cellmax({DATA.Trials(id).ytype});
    Expt.Stimvals.e3 = cellmax({DATA.Trials(id).ztype});
    Expt.Header.BlockStart = E.firsttrial;
    Expt.Header.BlockStartId = Expt.Trials(1).id;
    Expt.Header.rc = 0;
    Expt = FillStimvals(Expt);
    if isfield(DATA,'filename')
        Expt.Header.Name = DATA.filename;
    end
    Expt.Header = CopyFields(Expt.Header,DATA,'CreationDate');
        
function oExpts = MakeExpts(DATA, Expts)
    mintrials = 4;
    
    for j = 1:length(Expts)
        oExpts{j} = MakeOneExpt(DATA, Expts(j));
        if isfield(oExpts{j},'Trials') && length(oExpts{j}.Trials) > mintrials
            good(j) = 1;
        else
            good(j) = 0;
        end
    end
    oExpts = oExpts(good>0);
        
function s = cellmax(C)
    [a,b] = Counts(C);
    [~,c] = max(a);
    s = b{c};

function Expts = MakeAllExpts(DATA, varargin)

Expts = {};
psychonly = 0;
if ~isfield(DATA,'score')
    return;
end
[id, eid] = FindBlocks(DATA);
if psychonly
    sid = find(DATA.score <= 5);
else
    sid = find(DATA.score <= 5 | ismember(DATA.score,[51 50]));
end

%Check for unterminated blocks e.g. becuase of crash
for j = 1:length(id)    
    if j < length(id) && eid(j) > id(j+1)
        eid(j+1:end+1) = eid(j:end);
        eid(j) = id(j+1)-1;
    end
end
nx = 1;
for j = 1:length(id);
    Expt = [];
    tmp = DATA;
%this union seems to conflate all the blocks.  ? not sid intented   ?
%intersect? Changed to intersect Jan 2015. But tests in MakeExpt anyway
    ids = intersect([id(j):eid(j)],sid);
    ids = id(j):eid(j);
    tmp.score = DATA.score(ids);
    if isfield(DATA,'trialvals')
        tmp.trialvals = DATA.trialvals(ids);
    end
    tmp.x = DATA.x(ids);
    tmp.y = DATA.y(ids);
    if isfield(DATA,'ztype')
        tmp.z = DATA.z(ids);
        tmp.ztype = DATA.ztype(ids);
    end
    tmp.xtype = DATA.xtype(ids);
    tmp.ytype = DATA.ytype(ids);
    tmp.sign = DATA.sign(ids);
    tmp.times = DATA.times(ids);
    tmp.rwszs = DATA.rwszs(ids);
    tmp.rwsum = DATA.rwsum(ids);
    tmp.DURS = DATA.DURS(ids);
    tmp.saved = DATA.saved(ids);
    Expt = MakeExpt(tmp, varargin{:});
    f = fields(DATA.Block);
    n = (j * 2)-1;
    if isfield(Expt,'Trials') && length(Expt.Trials) > DATA.plot.mintrials
        goodexpt(j) = 1;
    else 
        goodexpt(j) = 0;
    end
    if DATA.score(eid(j)) == 27 && ~DATA.useallexpts
        goodexpt(j) = 0;
    end
    if goodexpt(j)
        for k = 1:length(f)
            bid = find(DATA.Blockid.(f{k}) <= id(j));
            if length(bid) & bid(end) <= length(DATA.Block.(f{k}))
                Expt.Stimvals.(f{k}) = DATA.Block.(f{k})(bid(end));
            else
                Expt.Stimvals.(f{k}) = NaN;
            end
        end
        Expts{nx} = Expt;
        if DATA.score(eid(j)) == 27
           Expts{nx}.Header.result = 19;
        else
           Expts{nx}.Header.result = 2;
        end
        nx = nx+1;
    end
end
if ~isempty(DATA.Comments)
    ts = 0;
    for j = 1:length(Expts)
        te = Expts{j}.Trials(end).End(end)+10000;
        cid = find([DATA.Comments.time] > ts & [DATA.Comments.time] <= te);
        Expts{j}.Comments.text = {DATA.Comments(cid).text};
        Expts{j}.Comments.times = [DATA.Comments(cid).time];
        ts = te;
    end
end

if DATA.verbose
    fprintf('%d Expts Made\n',nx);
end

function Expt = MakeExpt(DATA, varargin)

useall = 0;
skipblock = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'skipblock',5)
        j = j+1;
        skipblock = varargin{j};
    elseif strncmpi(varargin{j},'useall',5)
        useall = 1;
    end
    j = j+1;
end
if DATA.useallexpts
    sid = ismember(DATA.score, [1 2 51]);
    useall = 1;
else
    sid = ismember(DATA.score, [1 2]);
end
Expt = [];

if sum(sid) == 0 %no psych trials
    if useall
        sid = find(DATA.score ==51);
        if isempty(sid)
            return;
        end
    else
        return;
    end
end
[et, p] = GetType(DATA.xtype(sid));
[e2, p] = GetType(DATA.ytype(sid));
if isfield(DATA,'ztype')
[e3, p] = GetType(DATA.ztype(sid));
else
    e3 = 'e0';
end

Expt.Stimvals.et = et;
if ~isfield(Expt.Stimvals,'e2')
Expt.Stimvals.e2 = 'e0';
end

if ~isfield(Expt.Stimvals,'e3')
    if sum(strncmp(e3,{'se'},2))
        Expt.Stimvals.e3 = 'e0';
    else
        Expt.Stimvals.e3 = e3;
    end
end
chkf = {'ei' 'i2' 'i3'};
defaults = {0 0 0};
for j = 1:length(chkf)
if ~isfield(Expt.Stimvals,chkf{j})
    Expt.Stimvals.(chkf{j}) = defaults{j};
end
end

fn = fields(DATA.Stimvals);
for j = 1:length(fn)
    Expt.Stimvals.(fn{j}) = DATA.Stimvals.(fn{j});
end

Expt.Header.rc = 0;
Expt.Header.expname = DATA.filename;
if useall
    usetrials = [0 1 3 51];
else
    usetrials = [0 1 3];
end
id = find(ismember(DATA.score,usetrials));
if length(id) > 1 && (length(unique(DATA.y(id))) > 1 || DATA.Stimvals.ve > 5 || DATA.Stimvals.ve < 2)
    Expt.Stimvals.e2 = e2;
    blocks(1) = id(1);
else
    Expt = [];
    return;
end
nb = 1;
if isfield(DATA,'trialvals')
f = fields(DATA.trialvals);
else 
    f = [];
end

if DATA.verbose
    fprintf('Filling %d trials..',length(id));
end
if id(end) > length(DATA.y)
    fprintf('Y too short!!\n');
end
if id(end) > length(DATA.x)
    fprintf('X too short!!\n');
end
if id(end) > length(DATA.sign)
    fprintf('Sign too short!!\n');
end
Expt.Header.ntrials = length(id);
  for j = 1:length(id)
      if id(j) > 1 && DATA.score(id(j)-1) == 5
          blocks(nb) = id(j);
          nb = nb+1;
      end
      Expt.Trials(j).(et) = DATA.x(id(j));
      if ~strcmp(Expt.Stimvals.e2,'e0')
          Expt.Trials(j).(e2) = DATA.y(id(j));
      end
          
      Expt.Trials(j).rwdir = DATA.sign(id(j));
      if DATA.score(id(j)) == 1
          Expt.Trials(j).RespDir = DATA.sign(id(j));
      elseif DATA.score(id(j)) == 0
          Expt.Trials(j).RespDir = -DATA.sign(id(j));
      else
          Expt.Trials(j).RespDir = 0;
          Expt.Trials(j).rwdir = 0;
      end
      Expt.Trials(j).Start = DATA.times(id(j)).*10000;
      Expt.Trials(j).End = (DATA.times(id(j)).*10000)+20000;
      Expt.Trials(j).Trial = id(j);
      Expt.Trials(j).rw = DATA.rwszs(id(j));
      Expt.Trials(j).rwsum = DATA.rwsum(id(j));
      for k = 1:length(f)
          if length(DATA.trialvals(id(j)).(f{k}))
              Expt.Trials(j).(f{k}) = DATA.trialvals(id(j)).(f{k});
          end
      end
      if isfield(Expt.Trials,'TwoCylDisp') && isfield(Expt.Trials,'hx')
          Expt.Trials(j).dx = Expt.Trials(j).TwoCylDisp;
          Expt.Trials(j).bd = Expt.Trials(j).TwoCylDisp;
          if Expt.Trials(j).TwoCylDisp == -1011  %flip
              Expt.Trials(j).dx = abs(Expt.Trials(j).hx) .* DATA.sign(id(j));
              Expt.Trials(j).bd = -Expt.Trials(j).dx;
          elseif Expt.Trials(j).TwoCylDisp == 0  %flip
              Expt.Trials(j).bd = Expt.Trials(j).hx;
          end
          Expt.Trials(j).rd = Expt.Trials(j).dx - Expt.Trials(j).bd;
      end
      Expt.Trials(j).OptionCode = '+2a';
      Expt.Trials(j).Spikes = 0;
  end
  blocks(nb) = max([Expt.Trials.Trial]);
if DATA.verbose
    fprintf('..Done\n');
end
  if skipblock
      t = blocks(skipblock+1)-1;
      Expt.Header.BlockStart = blocks(skipblock+1:end);
      Expt.Trials = Expt.Trials(t:end);
  else
      Expt.Header.BlockStart = blocks;
  end
 id = find([Expt.Trials.rwdir] == 0);
 Expt.Header.Start = Expt.Trials(1).Start(1);
 Expt.Header.End = Expt.Trials(end).End(end);
 Expt.Header = CopyFields(Expt.Header,DATA,{'CreationDate'});
 %When reading an onlien file, don't know if smr files were split or not
 % treat them as if not, sp keep one creationdate, then references all
 % times to that.
 if 0
 id = find(DATA.score == 4 & DATA.times > 1000);
 if ~isempty(id)
     Expt.Header.CreationDate = DATA.times(id(1));
 end
 end
 
 function [type,frac] = GetType(x)

v = unique(x);
for j = 1:length(v)
    n(j) = length(strmatch(v{j},x));
end
[a,b] = max(n);
type = v{b};
frac = a./sum(n);

function dn = ctime2datenum(x)
         
    
    dn = x./(24 * 60 * 60);
%dont understand where 20/24 comes from. Someone doesn't use midnight    
    dn = dn + 719528 + (20/24); %add 1 Jan 1970
    