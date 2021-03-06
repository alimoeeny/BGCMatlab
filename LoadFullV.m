function FullV = LoadFullV(name, varargin)
%FullV = LoadFullV(name)
%loads a FullV file from disk, and performs any necessary scaling
%
%FullV = LoadFullV(name,'noconvert') leaves FullV.V as int16
%                  instead of converting to double
%FullV = LoadFullV(name,'forceall') forces inclusion of data marked
%           as not needed (e.g. when file was split

if iscell(name)
    for j = 1:length(name)
        res(j) = LoadFullV(name{j}, varargin{:});
    end
    return;
end

highpass = 1;
smoothw = 100;
keepsmooth = 0;
addtime = 0;
convert = 1;
silent = 0;
toint = 0;
loadmean = 0;
forceall = 0;
FullV = [];
meandata = [];
loadmeanfile = [];
if ~exist(name,'file')
    fprintf('%s Does not exist\n',name);
    return;
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'alldata',6)
        forceall = 1;
    elseif strncmpi(varargin{j},'addtime',6)
        addtime = 1;
    elseif strncmpi(varargin{j},'highpass',5)
        j = j+1;
        smoothw = varargin{j};
        keepsmooth = 1;
    elseif strncmpi(varargin{j},'checkint',6)
        toint = 3;
    elseif strncmpi(varargin{j},'converttoint',5)
        toint = 1;
    elseif strncmpi(varargin{j},'meanfile',5)
        j = j+1;
        loadmean = 1;
        meandata = varargin{j};
        if isfield(meandata,'name')
            loadmeanname = meandata.name;
        end
    elseif strncmpi(varargin{j},'saveint',5)
        toint = 2;
    elseif strncmpi(varargin{j},'silent',5)
        silent = 1;
    elseif strncmpi(varargin{j},'noconvert',5)
        convert = 0;
    elseif strncmpi(varargin{j},'nohighpass',5)
        highpass = 0;
    end
    j = j+1;
end


isint = 0;
ts = now;
try
    load(name);
catch ME
    CheckExceptions(ME);
    FullV.loadname = name;
    FullV.errmsg = {};
    FullV.loadtime = 0;
    FullV.size = 0;
    FullV.exptno = GetExptNumber(name);
    FullV.program = 'LoadFullV';
    FullV = AddError(FullV,'-show','Cannot Read %s',name);
    return;
end
FullV.loadname = name;
if ~isfield(FullV,'V')
    if ~isfield(FullV,'errmsg')
        FullV.errmsg = {};
    end
  FullV = AddError(FullV,'-show','Error loading %s no V data\n', name);
  return
end
if isfield(FullV,'errdata') && iscell(FullV.errdata)
    FullV.errdata = CellToStruct(FullV.errdata,'flat');
end

%If this fullv has been split, last sample records the max sample
%point before the split.  So all data is inside here for re-splitting
%but gets truncated here in routine load
if isfield(FullV,'lastsample') && forceall == 0
    b = FullV.lastsample;
    FullV.t = FullV.t(1:b);
    FullV.V = FullV.V(:,1:b);
    FullV.meanV = FullV.meanV(1:b);
end

if loadmean && isfield(FullV,'sumscale') %A Utah file with means removed
    mname = regexprep(name,'\.p[0-9]+FullV.mat','FullVmean.mat');
    if ~strcmp(loadmeanfile,mname)
        if exist(mname)
            load(mname);
            FullV.meandata.meanV = sumv;
            FullV.meandata.name = mname;
            if ~isfield(MeanV,'probes')
                FullV.meandata.probes = 1:96;
            else
                FullV.meandata.probes = MeanV.probes;
            end
        end
    end
end
if toint && isfloat(FullV.V)
    if isinteger(FullV.V) %alreay i sint
        isint = 1;
        fprintf('%s isalready ints\n',name);
        return;
    elseif toint == 3 %just checking
        fprintf('%s is double\n',name);
        return;
    end
    intscale(2) = double(intmax('int16')-5);
    intscale(1) = max(abs(FullV.V(:)));
     FullV.V = int16(round(double(FullV.V .* intscale(2)./intscale(1))));    
     FullV.intscale = intscale;
     if toint == 2
         fprintf('Saving %s as int\n',name);
         save(name,'FullV');
     end
    return;
elseif isfloat(FullV.V) %shouln't happen normally
    FullV.loadtype = 'double';
else
    FullV.loadtype = 'int';    
end

method = 3;  %change manually for testing
FullV.initialloadtime = mytoc(ts);
d = whos('FullV');
FullV.readrate = d.bytes./FullV.initialloadtime;
FullV.size = d.bytes./(1024.*1024);
if convert && isinteger(FullV.V)  && isfield(FullV,'intscale')
    if method == 1  %try to save memory
        np = size(FullV.V,1);
        for j = 1:np
            NewV(np+1-j,:) = double(FullV.V(np+1-j,:)) .* FullV.intscale(1)/FullV.intscale(2);
            FullV.V = FullV.V(1:np-j,:);
        end
        FullV.V = NewV;
    elseif method == 2
        FullV.V = double(FullV.V) .*FullV.intscale(1)/FullV.intscale(2);
    elseif method == 3
        FullV.V = double(FullV.V);
        np = size(FullV.V,1);
        for j = 1:np
            FullV.V(j,:) = FullV.V(j,:) .* FullV.intscale(1)/FullV.intscale(2);
        end
    else
        NewV = double(FullV.V) .*FullV.intscale(1)/FullV.intscale(2);
        FullV.V = NewV;
        clear NewV;
    end
end
FullV.loadtime = mytoc(ts);
if keepsmooth
    FullV.highpass = smoothw;
end
if isfield(FullV,'chspk') &&  isempty(strfind(name,'p96FullV')) && FullV.chspk(1) == 96    
    FullV.chspk = GetProbeFromName(name);
    mycprintf('red','\nSetting FullV.chspk to %d\n',FullV.chspk);
end
if highpass && isfield(FullV,'highpass') && ~isnan(FullV.highpass)
    ts = now;
    if FullV.highpass > 0
        smoothw = FullV.highpass;
    end
     sm = smooth(double(FullV.V),smoothw);
     if keepsmooth
         FullV.Vsmooth = sm;
     end
     FullV.V = int16(round(double(FullV.V) - sm));
     FullV.filtertime = mytoc(ts);
end
if addtime
    FullV.t = BuildFullVt(FullV);
end
if isfield(FullV,'firstblk') && FullV.firstblk > 1 && strfind(FullV.loadname,'aFullV');
    FullV.exptno = floor(FullV.exptno)+0.1;
end
if isfield(FullV,'errmsg')
    for j = 1:length(FullV.errmsg)
        if ~isfield(FullV,'errdata') || length(FullV.errdata) < j
            FullV.errdata(j).eid = FullV.exptno;            
        end
        if isfield(FullV.errdata,'time')
            dstr = sprintf('on %s',datestr(FullV.errdata(j).time));
        else
            dstr = '';
        end
        if ~silent
            cprintf('blue','FullV file said:%s%s\n',FullV.errmsg{j},dstr);
        end
        if ~isfield(FullV.errdata,'program') || isempty(FullV.errdata(j).program)
            FullV.errdata(j).program = 'LoadFullV';
        end
    end
end
for j = 1:length(FullV.blkstart)
    FullV.blkend(j) = FullV.blkstart(j)+ FullV.blklen(j).*FullV.samper;
end
