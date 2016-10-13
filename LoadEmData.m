
function [cExpt, details]= LoadEmData(cExpt, varargin)
%
%
%  [Expt, details] = LoadEmData(Expt) finds the matching .em matlab file, read in
% the eye position data, and adds these to each trial in Expt
%
% Adds a field EyeData to Each Trials
% Eyedata(:,1) = LH, 2 = RH, 3 = LV 4 = RV
%
%  LoadEmData(Expt,'reread')  forces re-reading of EM data from .em.mat
%  files
%  LoadEmData(Expt,'int')  makes EyeData ints (small for saving)
%  LoadEmData(Expt,'double')  makes EyeData doubles (needed for most
%  calculation)


global bgcfileprefix;

datatype = 'default';
details = [];
reread = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'int',3)
        if reread
        elseif isfield(cExpt,'Trials') 
            if isfield(cExpt.Trials,'EyeData')
                cExpt = FixEmTrials(cExpt, varargin{:});
                return;
            elseif strncmpi(varargin{j},'intconv',7) %only convert. do nothing if no data. For Spike2
                return;
            end
        end
        datatype = 'int';
    elseif strncmpi(varargin{j},'double',6)
        if reread
        elseif isfield(cExpt,'Trials') && isfield(cExpt.Trials,'EyeData')
            cExpt = FixEmTrials(cExpt, varargin{:});
            return
        elseif isfield(cExpt,'Trials') && isfield(cExpt.Trials,'Eyevals')
            cExpt = FixEmTrials(cExpt, varargin{:});
            return
        end
        datatype = 'double';
    elseif strncmpi(varargin{j},'reread',3)
        reread = 1;
    end
    j = j+1;
end

  if isfield(cExpt,'Expt') %An AllExpt
      cExpt.Expt = LoadEmData(cExpt.Expt,varargin{:});
      return;
  end
  if iscell(cExpt)
      for j = 1:length(cExpt)
          cExpt{j} = LoadEmData(cExpt{j},varargin{:});
      end
      return;
  end
      name = strrep(cExpt.Header.Name,'\','/');
  if isfield(cExpt.Header,'Spike2Version')
      [cExpt, details] = LoadSpike2EMData(cExpt, varargin{:});
      if strcmp(datatype,'double')
          cExpt = FixEmTrials(cExpt,'double');
      else
          cExpt = FixEmTrials(cExpt,'int');
      end
      return;
  else
      idx = regexp(name,'\.c[0-9]\.');
      emname = strrep(name,name(idx:idx+3),'.em.');
  end
  ok = 1;
  if ~exist(emname,'file')
      emname = [bgcfileprefix emname];
      if ~exist(emname,'file')
          fprintf('No EM file %s\n',emname);
          ok = 0;
      end
  end
  if ok
      load(emname);
      for j = 1:length(Expt.Trials)
          emlen(j) = length(Expt.Trials(j).Eyevals.lh);
      end
      emlen = min(emlen(find(emlen > 0)));
      cExpt.Header.emlen = emlen;
      for j = 1:length(Expt.Trials)
          if Expt.Trials(j).Start == cExpt.Trials(j).Start
              cExpt.Trials(j).EyeData = [];
              cExpt.Trials(j).EyeData(:,1) = Expt.Trials(j).Eyevals.lh(1:emlen);
              cExpt.Trials(j).EyeData(:,2) = Expt.Trials(j).Eyevals.rh(1:emlen);
              cExpt.Trials(j).EyeData(:,3) = Expt.Trials(j).Eyevals.lv(1:emlen);
              cExpt.Trials(j).EyeData(:,4) = Expt.Trials(j).Eyevals.rv(1:emlen);
              cExpt.Trials(j).goodem = 1;
          else
              cExpt.Trials(j).goodem = 0;
          end
      end
  end

  
function Expt = FixEmTrials(Expt, varargin)
j = 1;

toint = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'int',3)
        toint = 1;
    elseif strncmpi(varargin{j},'double',6)
        toint = 0;
    end
    j = j+1;
end

if toint
    if ~isfield(Expt.Header,'EMdataclass') || strcmp(Expt.Header.EMdataclass,'double')
        if ~isfield(Expt.Header,'emscale') && isfield(Expt.Trials,'EyeData')
        for j = 1:length(Expt.Trials)
            x = minmax(Expt.Trials(j).EyeData(~isnan(Expt.Trials(j).EyeData)));
            if length(x) == 2
                emrange(j,:) = x;
            end
        end
        emrange = minmax(emrange(:));
        emscale =  max(abs(emrange));
        Expt.Header.emscale(1) = emscale;
        maxint = double(intmax('int16')-5);
        Expt.Header.emscale(2) = maxint;        
        else
            emscale = Expt.Header.emscale(1);
            maxint = Expt.Header.emscale(2);
        end
        for j = 1:length(Expt.Trials)
            x = Expt.Trials(j).EyeData;
            x = round(x .*maxint./emscale);
            x(isnan(x)) = maxint+2;
            iEM{j} = int16(x);
            if isfield(Expt.Trials,'pupil')
                x = Expt.Trials(j).pupil;
                x = round(x .*maxint./emscale);
                x(isnan(x)) = maxint+2;
                Expt.Trials(j).pupil = int16(x);
            end
        end
        for j = 1:length(Expt.Trials)
            Expt.Trials(j).EyeData = iEM{j};
        end
        Expt.Header.EMdataclass = 'int';
    end
else
    if ~isfield(Expt.Trials,'EyeData')
        Expt = AddError(Expt,'-show','No EM Data in %s\n',GetName(Expt,'withsuffix'));
        return;
    end
    if ~isfield(Expt.Header,'EMdataclass') || strcmp(Expt.Header.EMdataclass,'int')
%if already double, return
        if ~isfield(Expt.Header,'EMdataclass') && ~isinteger(Expt.Trials(1).EyeData)
            Expt.Header.EMdataclass = 'double';
            return;
        end
        maxint = Expt.Header.emscale(2);
        emscale = Expt.Header.emscale(1);
        pad = NaN;
        if isfield(Expt.Trials,'EyeData')
            for j = 1:length(Expt.Trials)
                x = double(Expt.Trials(j).EyeData);
                x(x == maxint+2) = pad;
                dEM{j} = x .*emscale./maxint;
                if isfield(Expt.Trials,'pupil')
                    x = double(Expt.Trials(j).pupil);
                    x(x == maxint+2) = pad;
                    x = x .*emscale./maxint;
                    Expt.Trials(j).pupil = x;
                end
            end
            for j = 1:length(Expt.Trials)
                Expt.Trials(j).EyeData = dEM{j};
            end
        elseif isfield(Expt.Trials,'Eyevals')
            f = fields(Expt.Trials(1).Eyevals);
            for j = 1:length(Expt.Trials)
                for k = 1:length(f)
                    x = double(Expt.Trials(j).Eyevals.(f{k}));
                    Expt.Trials(j).Eyevals.(f{k}) = x .*emscale./maxint;
                end
            end
        end
        Expt.Header.EMdataclass = 'double';
    end
end
  
  
function [cExpt, details] = LoadSpike2EMData(cExpt, varargin)

selectile=5; 
name = strrep(cExpt.Header.Name,'\','/');
emname = strrep(name,'.mat','.em.mat');
if isunix && strncmp(name,'Z:',2)
    emname = strrep(emname,'Z:','/bgc');
elseif isunix && (strncmp(name,'/data',5) || strncmp(name,'/smr',4))
    emname = ['/b' emname];
end
fixid = 0;
reread = 0;
sacargs ={};
suffix = 0;
addsoft = 0;
verbose =1;
checknch = 1; 

usetrials = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'addsoft',5)
        addsoft = 1;
    elseif strncmpi(varargin{j},'name',3)
        j = j+1;
        emname = varargin{j};
    elseif strncmpi(varargin{j},'reread',3)
        reread = 1;
    elseif strncmpi(varargin{j},'datdir',6)
        j = j+1;
        datdir = varargin{j};
        emname = strrep(emname,fileparts(emname),datdir);        
    elseif strncmpi(varargin{j},'rmonoc',5)
        sacargs = {sacargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'lmonoc',5)
        sacargs = {sacargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'selectile',9)
        j = j+1;
        selectile = varargin{j};  %set length of record based on this percentile
    elseif strncmpi(varargin{j},'suffix',5)
        j = j+1;
        suffix = varargin{j};
    elseif strncmpi(varargin{j},'trials',5)
        j = j+1;
        usetrials = find(ismember([cExpt.Trials.Trial],varargin{j}));
    elseif strncmpi(varargin{j},'zeropad',5)
        sacargs = {sacargs{:} varargin{j}};
    end
    j = j+1;
end

if isfield(cExpt.Header,'loadname')
    [prefix, fname, fsuffix] = fileparts(cExpt.Header.loadname);
    [a,b,c] = fileparts(emname);
    newname = [prefix '/' b c];
    
    if isdir(prefix) && exist(newname)
        emname = newname;
    else
        newname = name2path(cExpt.Header.Name);
        if exist(newname)
            emname = newname;
        else
            fname = regexprep(fname,'\..*','');
            newname = [prefix '/' fname '.em.mat']; 
            if exist(newname)
                emname = newname;
            end
        end
    end
end


if isfield(cExpt.Header,'bysuffix') && cExpt.Header.bysuffix == 1    
    cExpt.Header.bysuffix =2;
    if ~isfield(cExpt.Header,'errs')
        cExpt.Header.errs = {};
    end
    tid = [cExpt.Trials.Trial];
    id = find(diff(tid) < 0);
    if ~isempty(id) 
        if length(id) == length(cExpt.Header.BlockStart) -1; %easy fix
            cExpt.Header.errs{end+1} = sprintf('Fixed Trials.Trial');
            cExpt = FixExpt(cExpt,'Trials');
        else
            cExpt.Header.errs{end+1} = sprintf('Trials.Trial is not consecutive - needs fixing');
        end
        cprintf('error','%s\n',cExpt.Header.errs{end});
    end
    if isfield(cExpt.Header,'BlockStart')
        for j = 1:length(cExpt.Header.BlockStart)
            if j == length(cExpt.Header.BlockStart)
                tid = cExpt.Header.BlockStart(j):cExpt.Trials(end).Trial;
            else
                tid = cExpt.Header.BlockStart(j):cExpt.Header.BlockStart(j+1)-1;
            end
            if isfield(cExpt.Header,'suffixes')
                suffix = cExpt.Header.suffixes(j);
            else
                suffix = cExpt.Header.Combineids(j)
            end
            [cExpt, details{j}] = LoadSpike2EMData(cExpt, varargin{:},...
                'trials',tid,'suffix',suffix);
            for k = 1:length(details{j}.errs)
                cExpt.Header.errs{end+1} = details{j}.errs{k};
            end                
        end
    else
        [cExpt, details] = LoadSpike2EMData(cExpt, varargin{:});
    end
    cExpt.Header.bysuffix =1;
    if isfield(cExpt.Trials,'EyeData') %success
        cExpt = FixEmTrials(cExpt,'double');
        cExpt = CheckEyeData(cExpt, sacargs{:});
        cExpt = CheckSaccades(cExpt,sacargs{:});
    end
    return;
end
details.errs = {};

if suffix
    emname = regexprep(emname,'\.[0-9]*\.em.mat',['.' num2str(suffix) '.em.mat']);
end

ok = 1;
if ~exist(emname,'file') & exist('bgcfileprefix','var')
    emname = [bgcfileprefix emname];
    if ~exist(emname,'file')
        fprintf('No EM file %s\n',emname);
        ok = 0;
    end
end

if ~exist(emname,'file') && isfield(cExpt.Header,'loadname')
    if isfield(cExpt.Header,'bysuffix') && cExpt.Header.bysuffix == 2
        a = fileparts(cExpt.Header.loadname);
        [b,c] = fileparts(cExpt.Header.loadname);
        c = regexprep(c,'\..*','');
        emname = [a '/' c '.' num2str(suffix) '.em.mat'];
    else
        emname = strrep(cExpt.Header.loadname,'.mat','.em.mat');
    end
end
details.emname = emname;
if isfield(cExpt.Trials,'EyeData') 
    if ~isempty(usetrials)
    elseif reread
        cExpt.Trials = rmfield(cExpt.Trials,'EyeData');
    else
        return;
    end
end


if isempty(usetrials)
    usetrials = 1:length(cExpt.Trials);
end

if ~exist(emname,'file')
    details.errs{end+1} = sprintf('Cannot find %s, needed for Trials %d - %d Id %d - %d\n',...
        emname,usetrials(1),usetrials(end),cExpt.Trials(usetrials(1)).id,cExpt.Trials(usetrials(end)).id);
    [cExpt.Trials(usetrials).goodem] = deal(0);
    [cExpt.Trials(usetrials).emstarti] = deal(NaN);
    cprintf('red','%s\n', details.errs{end})
    return;
end


if ok
    load(emname);
    if ~exist('Expt','var')
        if isempty(strfind('.em',emname))
            cprintf('red','Cannot Find EM file for %s\n',emname);
        else
            fprintf('red','%s Does not contatin EM Data\n',emname);
        end
     ok = 0;
    end
     
end
if ok
        
    
    if isfield(Expt.Header,'AllEmData')
        cExpt = expt.LoadRawEmData(cExpt, Expt,varargin{:});
        return;
    end

%emtpy start values mess up [Expt.Trials.Start] below
for j = 1:length(Expt.Trials)
    if isempty(Expt.Trials(j).Start)
        Expt.Trials(j).Start = NaN;
    end
    if ~isfield(Expt.Trials,'ftime') || isempty(Expt.Trials(j).ftime)
        Expt.Trials(j).ftime = NaN;
    end
    if isfield(Expt.Trials,'Eyevals') && ~isempty(Expt.Trials(j).Eyevals)
        V = Expt.Trials(j).Eyevals;
        if isinteger(V.rh) && isfield(Expt.Header,'emscale');
            f = fields(V);
            for k = 1:length(f)
                if isinteger(V.(f{k}))
                    V.(f{k}) = int2double(V.(f{k}),Expt.Header.emscale(1));
                end
            end
        end
        ls(1) = max(find(~isnan(V.rh)));
        ls(2) = max(find(~isnan(V.lh)));
        ls(3) = max(find(~isnan(V.rv)));
        ls(4) = max(find(~isnan(V.lv)));
        if isfield(V,'xv')
            ls(5) = max(find(~isnan(V.xh)));
            ls(6) = max(find(~isnan(V.xv)));
        end
        ls = min(ls);
        V.rh = V.rh(1:ls);
        V.lh = V.lh(1:ls);
        V.rv = V.rv(1:ls);
        V.lv = V.lv(1:ls);
        if isfield(V,'xv')
            V.xv = V.xv(1:ls);
            V.xh = V.xh(1:ls);
        end
        if isfield(V,'rpupil')
            V.pupil = V.rpupil;
        end
        Expt.Trials(j).Eyevals = V;
    end
end


    
    %
% calculation of emlen should also look at start - ftime, so that traces
% are alinged on stimulus on. 
for j = usetrials
    [diffs(j), trial(j)] = min(abs(cExpt.Trials(j).Start(1)-[Expt.Trials.Start]));
    if isfield(Expt.Trials,'Eyevals')
        if isempty(Expt.Trials(trial(j)).Eyevals)
            emlen(j) = 0;
        else
        emlen(j) = min([length(Expt.Trials(trial(j)).Eyevals.rh)...
           length(Expt.Trials(trial(j)).Eyevals.lh) ...
           length(Expt.Trials(trial(j)).Eyevals.rv) ...
           length(Expt.Trials(trial(j)).Eyevals.lv)]);
        end
    elseif isfield(Expt.Trials,'rh')
        emlen(j) = min([length(Expt.Trials(trial(j)).rh)...
            length(Expt.Trials(trial(j)).lh)...
            length(Expt.Trials(trial(j)).rv)...
            length(Expt.Trials(trial(j)).lv)]);
    else
        emlen(j) = 0;
    end
end
if isfield(Expt.Trials,'id') && length(unique([Expt.Trials.id]))
    fixid = 1;
end

%channels can differ in lenght by 1 sample, so take one less that shortest
%for lh
emlens = emlen;
if max(diffs./10) > 0.5 || verbose > 1
    fprintf('%s: Max Start mismatch %.1fms\n',emname,max(diffs./10));
end
id = find(emlen > 0 & diffs < 500);
emlen = floor(prctile(emlen(id),selectile)) -1;
if isempty(emlen) || sum(emlen > 0) == 0
    cprintf('red','No Useable Trials in %s\n',emname);
    return;
end
cExpt.Header.emlen = emlen;
%CRsamplerate (ibherited from BW was samples per tick)
cExpt.Header.CRrates = Expt.Header.CRrates;
if Expt.Header.CRrates(1) > 0
    cExpt.Header.CRsamplerate = 1/(10000 * Expt.Header.CRrates(1));
else
    cExpt.Header.CRsamplerate = 597.01./10000; %% rate for uStim files - where eye pos is in wrong channels
end

pre = [Expt.Trials.Start] - [Expt.Trials.ftime];
if isfield(Expt.Header,'preperiod')
prepts = Expt.Header.preperiod./cExpt.Header.CRsamplerate;
else
prepts = 100;
end
cExpt.Header.emtimes = ([1:emlen] .* 1./cExpt.Header.CRsamplerate) - mean(pre(find(~isnan(pre))));
minprepts = floor(min(pre) .* cExpt.Header.CRsamplerate);
if ~isfield(Expt.Trials,'Eyevals')
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).Eyevals.lh = Expt.Trials(j).lh;
        Expt.Trials(j).Eyevals.rh = Expt.Trials(j).rh;
        Expt.Trials(j).Eyevals.rv = Expt.Trials(j).rv;
        Expt.Trials(j).Eyevals.lv = Expt.Trials(j).lv;
    end
end

softoff = [0 0 0 0];
if ~isfield(Expt.Trials,'softoff')
    addsoft = 0;
end
neyech = 4; %default
if isfield(Expt.Trials,'Eyevals')
    details.allpretimes = pre;
    pre = [];
    for j = usetrials
        k = trial(j);
        pre(j) = cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime;
        if ~isempty(Expt.Trials(k).Eyevals)
            f = fields(Expt.Trials(k).Eyevals);
            nch(j) = sum(ismember(f,{'rh' 'rv' 'lh' 'lv' 'xh' 'xv'}));
        else
            nch(j) = 0;
        end
    end
    neyech = round(max(nch));
    details.pretimes = pre;    
    pretime = min(pre(pre > 0));
    minprepts = floor(pretime .* cExpt.Header.CRsamplerate);
    if minprepts < 1
        cprintf('red','EMExpt%d Missing preperiod\n',GetExptNumber(cExpt));
    end
    for j = usetrials
        k = trial(j);
        prept(j) = floor((cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime) .* cExpt.Header.CRsamplerate);
    end
    minprepts = floor(prctile(prept(usetrials),5))-1;
    if minprepts < 1
        cprintf('red','EMExpt%d Missing preperiod in %d trials\n',GetExptNumber(cExpt),sum(prept(usetrials) < 1));
    end
    for j = usetrials
        starti(j) = prept(j)-minprepts+1;
        if starti(j) <= 0
            lens(j) = emlens(j);
        else
            lens(j) = emlens(j) - starti(j);
        end
    end
    lens(lens<0) = 0;
    maxerr = max(emlen-lens(pre>0));
    emlen = floor(prctile(lens(usetrials),95))-1; 
%emlen is target length total     
    %prept(j)+emlen-minprepts = emlens(j);
%    emlen = emlen-maxerr;
    for j = usetrials
        k = trial(j);
        if addsoft && length(Expt.Trials(k).softoff) >= 4
            softoff = Expt.Trials(k).softoff;
            cExpt.Trials(k).softoff = softoff; 
        end
        if abs(Expt.Trials(k).Start(1) - cExpt.Trials(j).Start(1)) < 500 && starti(j) > -emlen
            EV = Expt.Trials(k).Eyevals;
            cExpt.Trials(j).EyeData = [];
            st = starti(j); %where to read from in this trials Eyevals
            dst = 1;  %Where to write to in aligned EyeData
            if st <= 0
                dst = 1-st;
                st = 1;
            end
            addpts = emlen-(lens(j)-starti(j));
            if isfield(Expt.Trials(k).Eyevals,'lh')
                if addpts > 0
                    Expt.Trials(k).Eyevals.lh(emlens(j)+1:emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.rh(emlens(j)+1:emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.lv(emlens(j)+1:emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.rv(emlens(j)+1:emlens(j)+addpts) = NaN;
                    if isfield(EV,'xh');
                        Expt.Trials(k).Eyevals.xv(emlens(j)+1:emlens(j)+addpts) = NaN;
                        Expt.Trials(k).Eyevals.xh(emlens(j)+1:emlens(j)+addpts) = NaN;
                    end
                    if isfield(EV,'pupil');
                        Expt.Trials(k).Eyevals.pupil(emlens(j)+1:emlens(j)+addpts) = NaN;
                    end
                end
                
                cExpt.Trials(j).EyeData(dst:emlen,1) = Expt.Trials(k).Eyevals.lh(st:emlen+st-dst)+softoff(1);
            cExpt.Trials(j).EyeData(dst:emlen,2) = Expt.Trials(k).Eyevals.rh(st:emlen+st-dst)+softoff(2);
            cExpt.Trials(j).EyeData(dst:emlen,3) = Expt.Trials(k).Eyevals.lv(st:emlen+st-dst)+softoff(3);
            cExpt.Trials(j).EyeData(dst:emlen,4) = Expt.Trials(k).Eyevals.rv(st:emlen+st-dst)+softoff(4);
            if isfield(EV,'xh');
                cExpt.Trials(j).EyeData(dst:emlen,5) = Expt.Trials(k).Eyevals.xh(st:emlen+st-dst);
                cExpt.Trials(j).EyeData(dst:emlen,6) = Expt.Trials(k).Eyevals.xv(st:emlen+st-dst);
            end
            if isfield(EV,'pupil');
                cExpt.Trials(j).pupil(dst:emlen) = Expt.Trials(k).Eyevals.pupil(st:emlen+st-dst);
            end
            if dst > prepts;
                cExpt.Trials(j).EyeData(1:dst-1,:) = NaN;
            end
            else
                fprintf('Missing data in Trial %d\n',k);
            end
            cExpt.Trials(j).goodem = 1;
            cExpt.Trials(j).ftime = Expt.Trials(k).ftime;
            cExpt.Trials(j).emstarti = 1+ minprepts - dst;
            if dst > prepts
                cprintf('red','Missing samples in Trial %d (id%d)\n',cExpt.Trials(j).Trial,cExpt.Trials(j).id);
            end
            if dst > 50
                fprintf('Missing %d samples for Trial %d\n',dst,j);
            end
            if fixid & Expt.Trials(k).id ~= cExpt.Trials(j).id;                
                Expt.Trials(k).id = cExpt.Trials(j).id;
                fixid = fixid+1;
            end
        else
            cExpt.Trials(j).EyeData = zeros(emlen,neyech);
            cExpt.Trials(j).goodem = 0;
            cExpt.Trials(j).emstarti = NaN;
        end
        cExpt.Trials(j).emdelay = Expt.Trials(k).Start - cExpt.Trials(j).Start(1);
        szs(j,:) = size(cExpt.Trials(j).EyeData);
        sds(j,:) = nanstd(cExpt.Trials(j).EyeData);
    end
else
    for j = usetrials
        k = trial(j);
        if abs(Expt.Trials(k).Start - cExpt.Trials(j).Start) < 300
            cExpt.Trials(j).EyeData = [];
            cExpt.Trials(j).EyeData(1:emlen,1) = Expt.Trials(k).lh(1:emlen)+softoff(1);
            cExpt.Trials(j).EyeData(1:emlen,2) = Expt.Trials(k).rh(1:emlen)+softoff(2);
            cExpt.Trials(j).EyeData(1:emlen,3) = Expt.Trials(k).lv(1:emlen)+softoff(3);
            cExpt.Trials(j).EyeData(1:emlen,4) = Expt.Trials(k).rv(1:emlen)+softoff(4);
            cExpt.Trials(j).ftime = Expt.Trials(k).ftime;
            cExpt.Trials(j).goodem = 1;
            sds(j,:) = std(cExpt.Trials(j).EyeData);
        else
            cExpt.Trials(j).EyeData = zeros(emlen,neyech);
            cExpt.Trials(j).goodem = 0;
        end
    end
end
cExpt.Header.emlen = emlen;
cExpt.Header.emtimes = ([1:emlen] .* 1./cExpt.Header.CRsamplerate) - pretime;
cExpt.Header.emChannels = {'lh' 'rh' 'lv' 'rv'}; 
cExpt.Header = CopyFields(cExpt.Header,Expt.Header,{'EyeOffsets' 'EyeGains'});
end

if checknch && neyech > 5
    msd = prctile(sds,90);
    bid = find(msd < 0.1)
    if min(msd(1:4)) > 0.5 && max(msd(5:6)) < 0.1
        cExpt = AddError(cExpt,'Removing flat extra eye channels 5 and 6');
        for j = 1:length(cExpt.Trials)
            cExpt.Trials(j).EyeData = cExpt.Trials(j).EyeData(:,1:4);
        end
    elseif ~isempty(bid)
        cExpt = AddError(cExpt,'Channels %s look to be flat',sprintf('%d ',bid));        
    end
end
if suffix == 0 && ok
    evar = mean(var(cat(3,cExpt.Trials(usetrials).EyeData)),3);
    if length(sacargs)
        cExpt = CheckSaccades(cExpt,sacargs{:});
    elseif (evar(1) > evar(2) * 50 && evar(3) > evar(4) * 10) || (evar(1) > evar(2) * 10 && evar(3) > evar(4) * 50)
        cExpt = CheckSaccades(cExpt,'lmonoc');
    elseif evar(2) > evar(1) * 50 && evar(4) > evar(3) * 50
        cExpt = CheckSaccades(cExpt,'rmonoc');
    else
        cExpt = CheckSaccades(cExpt);
    end
end

if fixid > 1 %% Ids were changed. Save .em file with Trial.id set correctly
    save(emname,'Expt');
end