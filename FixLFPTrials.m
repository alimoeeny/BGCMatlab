function Expt = FixLFPTrials(Expt, varargin)
%Expt = FixLFPTrials(Expt)
%Makes sure that All Trials.LFP fields have the same length, padding
%front and back with NaNs;
%Expt = FixLFPTrials(Expt,'double') converts int to doubles
%Expt = FixLFPTrials(Expt,'int') converts double back to int (do this
%before saving

MAXPROBES = 96;
pad = NaN;
needft = 0;
gotft = 0;
fixdouble = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'zeropad',6)
        pad = 0;
    elseif strncmpi(varargin{j},'double',5)
        fixdouble = 1;
        if strncmpi(varargin{j},'doubleconv',9)
            if isfield(Expt.Trials,'LFP') && ~isinteger(Expt.Trials(1).LFP)
                return;
            end
        end
    elseif strncmpi(varargin{j},'int',3)
        fixdouble = 2;
    elseif strncmpi(varargin{j},'ft',2)
        needft = 1;
    elseif strncmpi(varargin{j},'fixspike',6)
        Expt = FixLFPSpike(Expt,varargin{:});
        return;
    end
    j = j+1;
end

if isfield(Expt,'Spikes') && iscell(Expt.Spikes) && isfield(Expt,'Expt') %AllExpt
    Expt.Expt = FixLFPTrials(Expt.Expt,varargin{:});
    return;
end


if ~isfield(Expt.Trials,'LFP') %load faileed
    return;
end

try
    
if ~isfield(Expt.Trials,'lfptime') & isfield(Expt.Trials,'ftime') %spike2 raw data
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).lfptime = Expt.Trials(j).ftime;
    end    
end

if fixdouble ==2 %convert to int
    if isfield(Expt.Header,'LFPclass') && strcmp(Expt.Header.LFPclass,'int')
        return;
    end    
    
    if ~isfield(Expt.Header,'lfpscale')
        for j = 1:length(Expt.Trials)
            x = minmax(Expt.Trials(j).LFP(~isnan(Expt.Trials(j).LFP)));
            if length(x) == 2
                lfprange(j,:) = x;
            end
        end
        lfprange = minmax(lfprange(:));
        lfpscale =  max(abs(lfprange));
        Expt.Header.lfpscale = lfpscale;
        maxint = double(intmax('int16')-5);
        Expt.Header.lfpscale(2) = maxint;
    end
    if ~isfield(Expt.Header,'ftlfpscale') && isfield(Expt.Trials,'FTlfp')
        for j = 1:length(Expt.Trials)
            x = minmax(Expt.Trials(j).FTlfp(~isnan(Expt.Trials(j).FTlfp)));
            if length(x) == 2
                lfprange(j,:) = x;
            end
        end
        lfprange = minmax(lfprange(:));
        ftlfpscale =  max(abs(lfprange));
        Expt.Header.ftlfpscale = ftlfpscale;
        maxint = double(intmax('int16')-5);
        Expt.Header.ftlfpscale(2) = maxint;
    elseif isfield(Expt.Header,'ftlfpscale')
        ftlfpscale = Expt.Header.ftlfpscale;s
    end
    maxint = Expt.Header.lfpscale(2);

    for j = 1:length(Expt.Trials)
        x = Expt.Trials(j).LFP;
        x = round(x .*maxint./lfpscale);
        x(isnan(x)) = maxint+2;
        iLFP{j} = int16(x);
        if isfield(Expt.Trials,'FTlfp')
            x = Expt.Trials(j).FTlfp;
            x = round(x .*maxint./ftlfpscale);
            x(isnan(x)) = maxint+2;
            iFT{j} = int16(x);
        end
    end
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).LFP = iLFP{j};
        if isfield(Expt.Trials,'FTlfp')
            Expt.Trials(j).FTlfp = iFT{j};
        end
    end
    Expt.Header.LFPclass = 'int';
%if this is the first callto FixLFPTrials, proceed    
    if isfield(Expt.Header,'fixed') %if 
        return;
    end
elseif fixdouble && isfield(Expt.Header,'lfpscale')
    maxint = Expt.Header.lfpscale(2);
    for j = 1:length(Expt.Trials)
        x = Expt.Trials(j).LFP;
        dLFP{j} = (double(x).*Expt.Header.lfpscale(1)./maxint);
        dLFP{j}(x == maxint+2) = pad;
        if isfield(Expt.Trials,'FTlfp')
            x = Expt.Trials(j).FTlfp;
            dFT{j} = (double(x).*Expt.Header.ftlfpscale(1)./maxint);
            dFT{j}(x == maxint+2) = pad;
        end
    end
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).LFP = dLFP{j};
        if isfield(Expt.Trials,'FTlfp')
            Expt.Trials(j).FTlfp = dFT{j};
        end
    end
    Expt.Header.LFPclass = 'double';
%if this is the first callto FixLFPTrials, proceed    
    if isfield(Expt.Header,'fixed') %if 
        return;
    end
end

for j = 1:length(Expt.Trials)
    sizes(j,:) = size(Expt.Trials(j).LFP);
    if isempty(Expt.Trials(j).lfptime)
        starts(j) = 0;
    else
    starts(j) = Expt.Trials(j).Start(1) - Expt.Trials(j).lfptime;
    end
end

if max(sizes(:,2)) > MAXPROBES
    flip = 1;
    sizes = fliplr(sizes);
else
    flip = 0;
end
maxl = floor(prctile(sizes(:,1),90));
ts = median(starts);

%samplerate is the duration of a sample in units of 1/10 ms
if isfield(Expt.Header,'LFPsamplerate')
    samplerate = Expt.Header.LFPsamplerate .* 10000; %this is in sec
elseif isfield(Expt.Header,'CRsamplerate')
    samplerate = Expt.Header.CRsamplerate .* 1000;
    if samplerate < 1 %must have been in 1/10 of ms
        samplerate = samplerate * 10;
    end
else
    samplerate = 10;
end

%samplerate is samplerate in ms.

if isfield(Expt.Header,'preperiod')
    preperiod = Expt.Header.preperiod;
else
    preperiod = ts;
end
    

for j = 1:length(Expt.Trials)
    if flip
        Expt.Trials(j).LFP = Expt.Trials(j).LFP';
    end
    np(j) = floor((ts - starts(j)) ./(samplerate));
    if ~isempty(Expt.Trials(j).LFP)
        lfprange(j,:) = minmax(Expt.Trials(j).LFP(:));
    end

    if np(j) > 0
        pre = ones(np(j),sizes(j,2)).* pad;
        Expt.Trials(j).LFP = cat(1,pre, Expt.Trials(j).LFP);
        Expt.Trials(j).lfpvalid(1) = np(j)+1;
    elseif np(j) < 0
        t = -np(j);
        Expt.Trials(j).LFP =  Expt.Trials(j).LFP(1+t:end,:);
        Expt.Trials(j).lfpvalid(1) = 1;
    end
    npost = maxl-size(Expt.Trials(j).LFP,1);
    if npost > 0
        post = ones(npost, sizes(j,2)).* pad;
        Expt.Trials(j).LFP = cat(1,Expt.Trials(j).LFP,post);
        Expt.Trials(j).lfpvalid(2) = maxl-npost;
    elseif npost <= 0
        Expt.Trials(j).LFP = Expt.Trials(j).LFP(1:maxl,:);
        Expt.Trials(j).lfpvalid(2) = maxl;
    end
    nposts(j) = npost;
    newsizes(j,:) = size(Expt.Trials(j).LFP);
    if needft
       Expt.Trials(j).FTlfp = fft(Expt.Trials(j).LFP);
    elseif isfield(Expt.Trials,'FTlfp')
        if size(Expt.Trials(j).FTlfp,1) < maxl
            Expt.Trials(j).FTlfp(end+1:maxl,:) = NaN;
        end
    end
end
    Expt.Header.lfplen = min(newsizes(:,1));

Expt.Header.LFPtimes = (([1:Expt.Header.lfplen] .* samplerate) - (preperiod));
if needft
    Expt.Header.LFPfreq = [1:Expt.Header.lfplen] .* 1./(Expt.Header.lfplen .* Expt.Header.LFPsamplerate);
end
Expt.Header.fixed = 1;
if fixdouble ~= 1
    Expt = FixLFPTrials(Expt,varargin{:},'int');
else
    Expt.Header.LFPclass = 'double';
end
max(np);

catch ME
    fid = fopen('C:\Spike2\FixLFPErrs.txt','a')
    str = CheckExceptions(ME);
    fprintf(fid,'Error Fixing %d\n%s\n',Expt.Header.Name,str);
    fclose(fid);
end