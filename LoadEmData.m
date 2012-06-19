
function cExpt = LoadEmData(cExpt, varargin)
%
%
%  Expt = LoadEmData(Expt) finds the matching .em matlab file, read in
% the eye position data, and adds these to each trial in Expt
%
global bgcfileprefix;

  name = strrep(cExpt.Header.Name,'\','/');
  if isfield(cExpt.Header,'Spike2Version')
      cExpt = LoadSpike2EMData(cExpt, varargin{:});
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

  
  
function cExpt = LoadSpike2EMData(cExpt, varargin)

name = strrep(cExpt.Header.Name,'\','/');
emname = strrep(name,'.mat','.em.mat');
if isunix && strncmp(name,'Z:',2)
    emname = strrep(emname,'Z:','/bgc');
elseif isunix && (strncmp(name,'/data',5) || strncmp(name,'/smr',4))
    emname = ['/bgc' emname];
end
fixid = 0;
reread = 0;
sacargs ={};

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'name',3)
        j = j+1;
        emname = varargin{j};
    elseif strncmpi(varargin{j},'reread',3)
        reread = 1;
    elseif strncmpi(varargin{j},'rmonoc',5)
        sacargs = {sacargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'lmonoc',5)
        sacargs = {sacargs{:} varargin{j}};
    end
    j = j+1;
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
    emname = strrep(cExpt.Header.loadname,'.mat','.em.mat');
end
if isfield(cExpt.Trials,'EyeData') 
    if reread
    cExpt.Trials = rmfield(cExpt.Trials,'EyeData');
    else
        return;
    end
end

if ~exist(emname,'file')
    return;
end

if ok
    load(emname);

%emtpy start values mess up [Expt.Trials.Start] below
for j = 1:length(Expt.Trials)
    if isempty(Expt.Trials(j).Start)
        Expt.Trials(j).Start = NaN;
    end
    if ~isfield(Expt.Trials,'ftime') || isempty(Expt.Trials(j).ftime)
        Expt.Trials(j).ftime = NaN;
    end
end


    
    %
% calculation of emlen should also look at start - ftime, so that traces
% are alinged on stimulus on. 
for j = 1:length(cExpt.Trials)
    [diffs(j), trial(j)] = min(abs(cExpt.Trials(j).Start(1)-[Expt.Trials.Start]));
    if isfield(Expt.Trials,'Eyevals')
        emlen(j) = min([length(Expt.Trials(trial(j)).Eyevals.rh)...
           length(Expt.Trials(trial(j)).Eyevals.lh) ...
           length(Expt.Trials(trial(j)).Eyevals.rv) ...
           length(Expt.Trials(trial(j)).Eyevals.lv)]);
    else
        emlen(j) = min([length(Expt.Trials(trial(j)).rh)...
            length(Expt.Trials(trial(j)).lh)...
            length(Expt.Trials(trial(j)).rv)...
            length(Expt.Trials(trial(j)).lv)]);
    end
end
if isfield(Expt.Trials,'id') && length(unique([Expt.Trials.id]))
    fixid = 1;
end

%channels can differ in lenght by 1 sample, so take one less that shortest
%for lh
emlens = emlen;
fprintf('Max Start mismatch %.1fms\n',max(diffs./10));
emlen = min(emlen(find(emlen > 0 & diffs < 500))) -1;
if isempty(emlen)
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
cExpt.Header.emtimes = ([1:emlen] .* 1./cExpt.Header.CRsamplerate) - mean(pre(find(~isnan(pre))));
minprepts = floor(min(pre) .* cExpt.Header.CRsamplerate);

if isfield(Expt.Trials,'Eyevals')
for j = 1:length(cExpt.Trials)
    k = trial(j);
    pre(j) = cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime;
end
pretime = min(pre);
minprepts = floor(min(pre) .* cExpt.Header.CRsamplerate);
for j = 1:length(cExpt.Trials)
    k = trial(j);
    prepts = floor((cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime) .* cExpt.Header.CRsamplerate);
    starti = prepts-minprepts+1;
    lens(j) = emlens(j) - starti;
end

maxerr = max(emlen-lens);
emlen = emlen-maxerr;
for j = 1:length(cExpt.Trials)
    k = trial(j);
    if abs(Expt.Trials(k).Start(1) - cExpt.Trials(j).Start(1)) < 500
        prepts = floor((cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime) .* cExpt.Header.CRsamplerate);
        starti = prepts-minprepts+1;
        cExpt.Trials(j).EyeData = [];
        cExpt.Trials(j).EyeData(:,1) = Expt.Trials(k).Eyevals.lh(starti:emlen+starti-1);
        cExpt.Trials(j).EyeData(:,2) = Expt.Trials(k).Eyevals.rh(starti:emlen+starti-1);
        cExpt.Trials(j).EyeData(:,3) = Expt.Trials(k).Eyevals.lv(starti:emlen+starti-1);
        cExpt.Trials(j).EyeData(:,4) = Expt.Trials(k).Eyevals.rv(starti:emlen+starti-1);
        cExpt.Trials(j).goodem = 1;
        Expt.Trials(k).id = cExpt.Trials(j).id;
    else
        cExpt.Trials(j).EyeData(1:emlen,1) = zeros(emlen,1);
        cExpt.Trials(j).EyeData(1:emlen,2) = zeros(emlen,1);
        cExpt.Trials(j).EyeData(1:emlen,3) = zeros(emlen,1);
        cExpt.Trials(j).EyeData(1:emlen,4) = zeros(emlen,1);
        cExpt.Trials(j).goodem = 0;
    end
    cExpt.Trials(j).emdelay = Expt.Trials(k).Start - cExpt.Trials(j).Start(1);
end
else
for j = 1:length(cExpt.Trials)
    k = trial(j);
    if abs(Expt.Trials(k).Start - cExpt.Trials(j).Start) < 300
        cExpt.Trials(j).EyeData = [];
        cExpt.Trials(j).EyeData(:,1) = Expt.Trials(k).lh(1:emlen);
        cExpt.Trials(j).EyeData(:,2) = Expt.Trials(k).rh(1:emlen);
        cExpt.Trials(j).EyeData(:,3) = Expt.Trials(k).lv(1:emlen);
        cExpt.Trials(j).EyeData(:,4) = Expt.Trials(k).rv(1:emlen);
        cExpt.Trials(j).goodem = 1;
    else
        cExpt.Trials(j).EyeData(:,1) = zeros(emlen,1);
        cExpt.Trials(j).EyeData(:,2) = zeros(emlen,1);
        cExpt.Trials(j).EyeData(:,3) = zeros(emlen,1);
        cExpt.Trials(j).EyeData(:,4) = zeros(emlen,1);
        cExpt.Trials(j).goodem = 0;
    end
end
end
end
cExpt.Header.emlen = emlen;
cExpt.Header.emtimes = ([1:emlen] .* 1./cExpt.Header.CRsamplerate) - pretime;

evar = mean(var(cat(3,cExpt.Trials.EyeData)),3);
if length(sacargs)
    cExpt = CheckSaccades(cExpt,sacargs{:});
elseif evar(1) > evar(2) * 50 && evar(3) > evar(4) * 10
    cExpt = CheckSaccades(cExpt,'lmonoc');
elseif evar(2) > evar(1) * 50 && evar(4) > evar(3) * 50
    cExpt = CheckSaccades(cExpt,'rmonoc');
else
    cExpt = CheckSaccades(cExpt);
end


if fixid %% Save .em file with Trial.id set correctly
    save(emname,'Expt');
end