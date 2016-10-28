function [Spikes, spkfile] = ReadSpikeFile(spkfile, varargin)
%Spikes = ReadSpikeFile(spkfile, varargin)
%    ....,'allprobes')  loads in other probes from xspk file
%    ....,'double')  converts int to double 
% if S is a Cluster struture, loads matching file...

loadallspikes = 0;
converttodouble = 0;
j = 1;
while j  <= length(varargin)
    if strncmpi(varargin{j},'allprobes',6)
        loadallspikes = 1;
    elseif strncmpi(varargin{j},'double',6)
        converttodouble = 1;
    end
    j = j+1;
end
spkpath = ['Spikes/' GetMonkeyName(spkfile)];

if iscell(spkfile)
    for j = 1:length(spkfile)
        [Spikes{j}, spkfiles{j}] = ReadSpikeFile(spkfile{j},varargin{:});
    end
    spkfile = spkfiles;
    return;
end
if isstruct(spkfile)
    S = spkfile;
    spkfile = BuildFileName(S,'spkfile');
else
    S = [];
end

if ~exist(spkfile)
    spkfile = strrep(spkfile,['Spikes/' GetMonkeyName(spkfile)],'Spikes/');
    if ~exist(spkfile)
        fprintf('Cannot read %s\n',spkfile);
        if ~isempty(strfind(spkfile,'Ecker'))
            fprintf('Use / in pathname  to read Ecker AutoSpikes\n');
        end
        Spikes = [];
        return;
    end
end
ts = now;
a = load(spkfile);
if isfield(a,'Spikes')
    Spikes = a.Spikes;
else
    f = fields(a);
    if length(f) ==1
        Spikes = a.(f{1});
    end    
end
%Spikes.loadname = spkfile;
if size(Spikes.values,2) > 100
    Spikes.values = Spikes.values';
end
Spikes.Header.error = 0;
Spikes.Header.loadname = spkfile;
    if isfield(S,'clst') && size(Spikes.values,1) ~= length(S.clst)
        Spikes = AddError(Spikes,'-show','Spike events (%d) do not match Cluster (%d)\n',size(Spikes.values,1),length(S.clst));
        Spikes.Header.error = 1;
    end
    Spikes.times = Spikes.times .* 10000;
    Spikes.times = reshape(Spikes.times,length(Spikes.times),1);
    if size(Spikes.codes,2) == 1
        Spikes.codes(:,2) = Spikes.codes(:,1);
    end
    Spikes.Header.bysuffix = 0;
    if ~isfield(Spikes.Header,'matfile')
%WE don't care about this any more it seems to me.
%When does Spikes.Header.bysuffix get used
%        fprintf('SpkFile %s missing Header.matfile\n',spkfile);
        Spikes.Header.matfile = '';
    end
    if regexp(Spikes.Header.matfile,'\.[0-9]*\.mat')
        Spikes.Header.bysuffix = 1;
    end
        
    if ~isfield(Spikes,'maxint')
        Spikes.maxint = double(intmax('int16')-5);
        if strmatch(class(Spikes.values),'int16')
            if ~isfield(Spikes,'maxv')
                fprintf('Integer Spikes without maxv!!\n');
                Spikes.maxv = double(max(abs(Spikes.values(:)))./5);
            end
            Spikes.VRange = double([min(Spikes.values(:))  max(Spikes.values(:))]).*Spikes.maxv./Spikes.maxint;
        else
            Spikes.maxint = 1;
            Spikes.maxv = max(abs(Spikes.values(:)));
            Spikes.VRange = [min(Spikes.values(:))  max(Spikes.values(:))];
        end
    else
        Spikes.VRange = double([min(Spikes.values(:))  max(Spikes.values(:))]).*Spikes.maxv./Spikes.maxint;
    end
    Spikes.Vscale = Spikes.maxv./Spikes.maxint;
    if isfield(Spikes,'VRange') && isempty(Spikes.VRange)
        
    end
if converttodouble && isinteger(Spikes.values)
    Spikes.values = double(Spikes.values) .* Spikes.Vscale;
end
if loadallspikes
    xname = regexprep(spkfile,'.p([0-9])*t','.p$1xt');
    if exist(xname)
        X = load(xname);
        Spikes.xmaxv = X.Spikes.maxv;
        Spikes.xVscale = X.Spikes.maxv./X.Spikes.maxint;
        if ndims(X.Spikes.values) == 2
            Spikes.xvalues(1,:,:) = X.Spikes.values;
        else
            Spikes.xVranges(:,1) = min(min(X.Spikes.values,[],3),[],2);
            Spikes.xVranges(:,2) = max(max(X.Spikes.values,[],3),[],2);
            Spikes.xVranges = Spikes.xVranges .* Spikes.xVscale;
            Spikes.xvalues = X.Spikes.values;
        end
        if ~isfield(X.Spikes,'xVrange')
            Spikes.xVrange = double([min(X.Spikes.values(:)) max(X.Spikes.values(:))]) .* X.Spikes.maxv./X.Spikes.maxint;
        else
            Spikes.xVrange = X.Spikes.xVrange;
        end
        Spikes.xchans = X.Spikes.chspk;
        if isfield(X.Spikes,'TriggerV')
            Spikes.TriggerV = X.Spikes.TriggerV;
        end
        if converttodouble && isinteger(Spikes.xvalues)
            Spikes.xvalues = double(Spikes.xvalues) .* Spikes.xVscale;
        end

    end
end
if ~isfield(Spikes,'exptno')
    Spikes.Header.exptno = GetExptNumber(spkfile);
end
Spikes.Header.loaddur = mytoc(ts);
Spikes.Header.loadtime = ts;
if isfield(S,'probe') && ~isfield(Spikes,'probe')
    Spikes.probe = S.probe;
end
