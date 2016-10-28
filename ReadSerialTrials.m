function [E, txt, details] = ReadSerialTrials(name,varargin)
%[Trials, txt] = ReadSerialTrials(name) Read serial file for AplaySpkfile
%Used for reading .bnc files
% convert a serial output file made by binoc in a Trials struct
% like ReadSerialFile but does not make expt strutures, and includes all
% fields
% ...,'savemat')  saves out a bnc.mat file
method = 2;
maxchar = 100;
verbose = 0;
details = [];
BlockEnds = [];
BlockStarts = [];
nexp = 1;
%longfields catches any codes > 2, so that rest can be grabbed blindly

longfields = {'expname'  'imve' 'exvals'  ...
     'Electrode' 'dxseq' 'ceseq'  'prev'};

longnfields = {'puA' 'puF'  'USp' 'USd' 'USf' 'imve'  'bjv' ...
    'nbars' 'ijump' 'mixac' 'imi'  'psyv' 'seof' 'serange' 'Trw' 'prev'};

specialfields = {'Reopened' 'EndExpt' 'EndStim' 'mt' 'BadFix' 'RW' 'RG' 'RF' 'RB' 'RL' 'Stimulus'};

ignorefields = {'by binoc'  'NewConnection' 'bpos' 'sb-' 'sb+' 'Experimenter' 'Hemifield' 'StartDepth' ...
    'vebinoclean' 'Sa:' 'HemiSphere' 'CLOOP' 'Hemisphere' 'Adapter' 'uf' 'O ' 'Log:' 'Remaining' 'cm=' 'stim:'};

charfields = {'et' 'st' 'e2' 'e3' 'Bs' 'lo' 'vve' 'cx' 'VisualArea' 'expvars'}; 
savemat = 0;


j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'newmethod',6)
        method = 2;
    elseif strncmpi(varargin{j},'savemat',6)
        savemat = 1;
    end
    j = j+1;
end

if ischar(name)
    txt = scanlines(name);
elseif iscellstr(name);
    txt = name;
end

t = 0;

Stimulus.Start = 0;
Stimulus.cellfield.mtFi = 0;
Stimuli = {};
tstart = now;
if isempty(txt)
    if ischar(name)
        cprintf('red','Cannot read %s\n',name);
    end
    return;
end
ignoreid = [];
for j = 1:length(ignorefields)
    id = find(strncmp(ignorefields{j},txt,length(ignorefields{j})));
    ignoreid = cat(1,ignoreid, id);
end
ignoreid = unique(ignoreid);
ignoreid(end+1) = length(txt)+1;
nextignore = 1;
txt = txt(setdiff(1:length(txt),ignoreid));
specialid = [];
for j = 1:length(specialfields)
    id = find(strncmp(specialfields{j},txt,length(specialfields{j})));
    specialid = cat(1,specialid, id);
end
nextspecial = 1;

charid = [];
for j = 1:length(charfields)
    id = find(strncmp(charfields{j},txt,length(charfields{j})));
    charid = cat(1,charid, id);
end
nextchar = 1;

longid = [];
for j = 1:length(longfields)
    id = find(strncmp(longfields{j},txt,length(longfields{j})));
    longid = cat(1,longid, id);
end
longid = unique(longid);

longidn = [];
for j = 1:length(longnfields)
    id = find(strncmp(longnfields{j},txt,length(longnfields{j})));
    longidn = cat(1,longidn, id);
end
longidn = unique(longidn);

allid = unique(cat(1,longid, specialid, charid,longidn));
allid(end+1) = length(txt)+1;

nextall = 1;
nextlong = 1;
nextignore = 1;
E.filename = name;
E.Trials = [];
fprintf('Loading text took %.2f\n',mytoc(tstart));
tstart = now;
E.cellfields = [];
E.vectorfields = [];
lastid = 0;
Stimulus.idoffset = 0;
lastprep = 0;
instim = 0;
nchar = 0;
lastresult = 0;
codes = zeros(1,length(txt));
codes(longid) = 'L';
codes(charid) = 'C';
codes(specialid) = 'S';
codes(longidn) = 'N';
%    hashid = find(strncmp('#',txt,1));
%    codes(hashid) = 'H';

Stimulus.mtFi = [];
Stimulus.id = 0;
opendates = [];

% for each trial outcome line ('RG' 'RB' 'RF' 'RW')
%for nper trials, dont get RF after each stim, so use 'EndStim'
%keep reading until the next #prep line (Start of next stim)
%Can sometimes get #prep twice with no RG in between. When this 
%happens it all belongs to the next R[GWB] line
for j = 1:length(txt);
    s = txt{j};
    code = codes(j);
    if length(s) < 3
    elseif s(1) == '#'
        if sum(strncmp(s,{'#Prep' '#Resetting'},5))  && instim %start of info describing next stimulus or end of expt block
            instim = 0;
%Some care is required for first trial. But not totally sure why we don't allow t == 0            
            if t == 0 && (j-lastprep) > 6
                t = 1;
            end
            if t > 0 && (j - lastprep) > 6 %if made from online, this line might be repeated
                if Stimulus.id < lastid
                    Stimulus.idoffset = lastid + Stimulus.idoffset;
                elseif Stimulus.id == lastid %repeat should not happen
%repeats can happen when a bad fix is at trial end, so that RB is followed by RF very quicky.  Here don't have to worry
%Trial with RF has correct id.
                    if verbose || Stimulus.result > 0
                        if lastresult > 0 || verbose
                            fprintf('Repeated id %d on line %d. Last prep line %d result %d (last%d)\n',lastid,j,lastprep,Stimulus.result,lastresult);
                        end
                    end
                end
                lastid = Stimulus.id;
                if method == 2
                    Stimuli{t} = rmfield(Stimulus,'cellfield');
                    if isfield(Stimulus,'Nf') && ~isempty(Stimulus.Nf) && length(Stimulus.mtFi) > Stimulus.Nf
                        Stimulus.Nf = [];
                    end
                    Stimulus.mtFn = [];
                    Stimulus.mtFl = [];
                    Stimulus.mtFi = [];
                    Stimulus.mtrS = [];
                    Stimulus.Nf = [];
                    if isfield(Stimulus,'ic')
                        Stimulus.cL = Stimulus.co+Stimulus.ic;
                        Stimulus.cR = Stimulus.co-Stimulus.ic;
                    else
                        Stimulus.cL = 1;
                        Stimulus.cR = 1;
                    end
                    f = fields(Stimulus.cellfield);
                    for k = 1:length(f)
                        if Stimulus.cellfield.(f{k}) > 1
                            Stimulus.(f{k}) = [];
                        end                            
                    end
                else
                f = setdiff(fields(Stimulus),'cellfield');
                for k = 1:length(f)
                    E.Trials(t).(f{k}) = Stimulus.(f{k});
                    if length(Stimulus.(f{k})) > 5 && ~ischar(Stimulus.(f{k}))
                        Stimulus.(f{k})  = []; %don't allow values to carry over
                    elseif sum(strncmp(f{k},{'mtFn' 'mtFl' 'mtFi' 'mtrS'},4))
                        Stimulus.(f{k})  = []; %don't allow values to carry over                            
                    end
                end
                end
                if isfield(Stimulus,'ce') && length(Stimulus.ce) > 20
                    Stimulus.ce = 1; %reset to ensure that ce seq doesn't persist
                end
                t = t+1;
                lastresult = Stimulus.result;
            end
            if t == 0
                t = 1;
            end
            if isfield(Stimulus,'mtFn') && ~isempty(Stimulus.mtFn)
%                diff(Stimulus.Fi);
            end
            lastprep = j;
            [a, b, c] = sscanf(s,'#Prep%d %d %f');
        elseif strncmp(s,'#du',3)  %start of info describing next stimulus
            a = sscanf(s,'#du%f(%f:%f');
            Stimulus.duration = a(1);
            if length(a) > 1
                Stimulus.nframes = a(2);
            end
            if length(a) > 2
                Stimulus.nomdur = a(3);
            end
        elseif strncmp(s,'#EXPTOVER',7) 
            BlockEnds(nexp) = t-1;
        end
    elseif j == allid(nextall)
        nextall = nextall+1;
        if code == 'L'
            id = find(strncmp(s,longfields,3));
            f = longfields{id(1)};
            Stimulus.(f) = s(length(f)+1:end);
            Stimuli = CheckStimuli(Stimuli,f,method);
            if isfield(Stimulus,'EndExpt')
                f = longfields{id(1)};               
            end
        elseif code == 'N'
            id = find(strncmp(s,longnfields,3));
            f = longnfields{id(1)};
            id = length(f)+1;
            if s(id) == '='
                id = id+1;
            end
            Stimulus.(f) = sscanf(s(id:end),'%f');
            Stimuli = CheckStimuli(Stimuli,f,method);
        elseif code == 'S'
            if strncmp(s,'EndStim',6)
                Stimulus.result = 2;
                instim = 1;                
            elseif strncmp(s,'BadFix',6)
                Stimulus.result = 0;
                instim = 1;                
            elseif strncmp(s,'Stimulus',6) %Start of New Expt
                nexp = nexp+1;
                BlockStarts(nexp) = t;
                BlockStartid(nexp) = Stimulus.id;
                Stimulus.expvars = ''; %Dont carry over
            elseif strncmp(s,'EndExpt',6)                
                nexp = nexp+1;
                BlockEnds(nexp) = t;
                BlockEndid(nexp) = Stimulus.id;
            elseif strncmp(s,'mtFl=',5)
                Stimulus.mtFl = sscanf(s(6:end),'%f');
                E.cellfields.mtFl = 1;
            elseif strncmp(s,'mtFn=',5) && length(s) > 5
                Stimulus.mtFn = sscanf(s(6:end),'%f');
                E.cellfields.mtFn = 1;
            elseif strncmp(s,'mtFi=',5)
                Stimulus.mtFi = sscanf(s(6:end),'%f');
                E.cellfields.mtFi = 1;
            elseif strncmp(s,'mtrS=',5)
                Stimulus.mtrS = sscanf(s(6:end),'%f');
                E.cellfields.mtrS = 1;
            elseif strncmp(s,'Reopened',6)
                if strncmp(s,'Reopened (bnc)',13) %Start of New Expt
                    s = s(20:end);
                    if isfield(Stimulus,'expvars')
                        Stimulus.expvars = '';
                    end
                elseif strncmp(s,'Reopened (uf)',13)
                    s = s(19:end);
                else
                    s = s(14:end);
                end
                opendate = my.datenum(s);
                if ~isnan(opendate)
                    opendates(end+1) = opendate;
                end
            elseif sum(strncmp(s,{'RG' 'RW' 'RB' 'RF'},2))
                if strncmp(s,'RG',2)
                    Stimulus.result = 1;
                elseif strncmp(s,'RW',2)
                    Stimulus.result = 1;
                elseif strncmp(s,'RB',2)
                    Stimulus.result = 0;
                elseif strncmp(s,'RF',2)
                    Stimulus.result = 1;
                end
                instim = 1;
                id = strfind(s,'id=');
                if ~isempty(id)
                    stimid = sscanf(s(id(1):end),'id=%d');
                    if stimid > Stimulus.id
                        if lastresult > 0
%in 4per stim 1 might finish, then stim 2 is prepared, but badfix happens before
%stim2 is presented. So id never gets sent

                            if Stimulus.result > 0
                            fprintf('id %d at end of trial, but %d at start\n',stimid,Stimulus.id);
                            end
                        end
                        Stimulus.id = stimid;
                    end
                end
                id = strfind(s,'bt=');
                if ~isempty(id)
                    ts = sscanf(s(id(1):end),'bt=%f');
                    if ts > 0
                        Stimulus.Start = ts;
                        Stimulus.TrialTime = opendate + ts./(60 * 60 * 24);
                    end
                end
                
            end
        elseif code == 'C'
            id = strfind(s,'=');
            if ~isempty(id)
                f = s(1:id(1)-1);
                nc = id(1)+1;
                Stimulus.(f) = s(nc:end);
            else
            end
        end
   elseif strfind(s,':')
        if ~isempty(strfind(s,'Reopened'))
        elseif ~isempty(strfind(s,'Resetting'))
        elseif ~isempty(strfind(s,'Run ended'))
        elseif strncmp(s,'Codes:',5)    
        elseif strncmp(s,'Outcodes:',5)    
        elseif strncmp(s,'epos:',5)    
        elseif strncmp(s,'Unrecognized',9)    
        else
        id = strfind(s,':');
        f = s(1:id(1)-1);
        if isvarname(f)  %If not a valid field name his is not a real sequence parameter
            E.cellfields.(f) = 1;
            a = sscanf(s(id(1)+1:end),'%f');
            Stimuli = CheckStimuli(Stimuli,f, method);
            Stimulus.(f) = a;
            Stimulus.cellfield.(f) = length(a);
        else
            fprintf('Ignoring %s\n',s);
        end
        end
%    elseif strfind(s,'=')
%        id = strfind(s,'=');
    elseif strfind(s,'=')
            id = strfind(s,'=');
            f = s(1:id(1)-1);
            nc = id(1)+1;
            x = sscanf(s(nc:end),'%f');
            if code == 'C'
                Stimulus.(f) = s(nc:end);                
            elseif ~isempty(x)
                Stimuli = CheckStimuli(Stimuli,f, method);
                try
                Stimulus.(f) = x;
                end
            end
    else
        if sum(strncmp(s,{'imx' 'imy'},3))
            f = s(1:3);
            nc = 4;
        else
            f = s(1:2);
            nc = 3;
        end
        if strcmp(f,'ce') && isfield(Stimulus,'ce') && length(Stimulus.ce) > 20
            Stimulus.ceval = sscanf(s(3:end),'%f');
        else
            if isvarname(f)
                try %in case of invalid fields
                    Stimuli = CheckStimuli(Stimuli,f, method);
                    Stimulus.(f) = sscanf(s(nc:end),'%f');
                end
            end
        end
        if isfield(Stimulus,f)
        if length(Stimulus.(f)) > 1
            if sum(strcmp(f,{'sq' 'vs' 'fp'})) %Only want first value
                Stimulus.(f) = Stimulus.(f)(1);
                Stimuli = CheckStimuli(Stimuli,f, method);
            else
                E.vectorfields.(f) = 1;
            end
        end
        if sum(strncmp(s,{'mtFi' 'mtFn'},4))
%                E.vectorfields.(f) = 1; %cell field in Trials.xx for idx
        end
            
        if isempty(Stimulus.(f)) %usually means sscanf above was char, not number
            Stimulus.(f) = s(nc:end);
            Stimuli = CheckStimuli(Stimuli,f, method);
            if code ~= 'C'  && nchar < maxchar
                fprintf('%s\n',s);
                nchar = nchar+1;
            end
        end
        end
    end
end
details.Stimuli = Stimuli;
details.BlockEnd = BlockEnds;
details.BlockStarts = BlockStarts;
E.Header.BlockStarts = BlockStarts;
details.opendates = opendates;
E.readdate = now;
if method == 2
    mytoc(tstart);
    ts = now;
    CheckStimuli(Stimuli,'xo',method);
    E.Trials = Stim2Trials(Stimuli);
    E = expt.Condense(E);
    E = expt.Condense(E,'fields');
end

if savemat && ischar(name)
    matfile = [name '.mat'];
    BinocExpt = E;
    fprintf('Saving %s\n',name);
    save(matfile,'BinocExpt');
end


function S = CheckStimuli(S, f, method)

if method < 2
    return;
end
if 0 %turn on to check freq of calls. 
persistent codes;


if ~isfield(codes,f)
    codes.(f) = 1;
else
    codes.(f) = codes.(f) + 1;
end
end


if isempty(S)
    return;
end

if ~isfield(S{end},f) && isvarname(f)
    fprintf('Adding Field %s\n',f);
    for j = 1:length(S)
        S{j}.(f) = 0;
    end
end

function Trials = Stim2Trials(S)

Stim = S{end};
allf = fields(Stim);
[f, ia] = setdiff(allf,fields(S{1}));
if ~isempty(f)
    for j = 1:length(S)
        [f, ia] = setdiff(allf,fields(S{j}));
        f = allf(sort(ia));
        for k = 1:length(f)
            S{j}.(f{k}) = 0;
        end
    end
end
for j = length(S):-1:1
    try
    Trials(j) = S{j};
    catch
        for k = 1:length(allf)
            Trials(j).(allf{k}) = S{j}.(allf{k});
        end
    end
end
