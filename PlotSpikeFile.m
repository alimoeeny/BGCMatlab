function result = PlotSpikeFile(spkfile, varargin)
% PlotSpikeFile(file, varargin)
% Plots the spike waveform data from a single file on disk
% useful for looking at spikes made by autocutting in AllVPcs
% PlotSpikeFile(file, 'quickspks')
Spks = [];
AllSpikes = {};
preperiod = 500;
postperiod = 500;
result = [];
Expt = [];
C = [];
excludetrials = [];
probe =1;
files{1} = spkfile;
quickspks = 0;
distancemeasure = '2Gauss';
parentfigure = [];
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j}, 'values')
            Spks = varargin{j};
        elseif isfield(varargin{j},'Trials')
            Expt = varargin{j};
        elseif isfield(varargin{j},'probe')
            C = varargin{j};
            probe = C.probe;
        end
    elseif strncmp(varargin{j},'distancemeasure',12)
        j = j+1;
        distancemeasure = varargin{j};
    elseif strncmp(varargin{j},'exclude',6)
        j = j+1;
        excludetrials = varargin{j};
    elseif strncmp(varargin{j},'parentfigure',6)
        j = j+1;
        parentfigure = varargin{j};
    elseif strncmp(varargin{j},'quickspks',6)
        quickspks = 100;
    elseif strncmp(varargin{j},'spkfile',6)
        j = j+1;
        files{length(files)+1} = varargin{j};
    elseif strncmp(varargin{j},'Trial',6)
%select trial. If used with 'quickspks', shows block with this trial
        j = j+1;
        selecttrial = varargin{j};
    end
    j = j+1;
end


if isfield(spkfile,'Spikes') %passed the Spks file, not name
    Spikes = spkfile.Spikes;
    Clusters = spkfile.Clusters;

    spkfile = Spikes.Header.loadname;
    [probe, ex, xs ] = Name2Probe(spkfile);
    probes = probe;
else
    [probes(1), ex, xs ] = Name2Probe(spkfile);
    if length(files) > 1
        for j = 1:length(files)
            probes(j) = Name2Probe(files{j});
            AllSpikes{probes(j)} = ReadSpikeFile(files{j});
        end
        probe = probes(1);
    else
        probe = probes(1);
        Spikes = ReadSpikeFile(spkfile);
        result.Spikes = Spikes;
        if isfield(Spikes,'probe')
            probe = Spikes.probe;
        end
    end
    cfile = strrep(spkfile,'/Spikes','');
    cfile = sprintf('%s/Expt%d%sClusterTimesDetails.mat',fileparts(cfile),ex,xs);
    afile = sprintf('%s/Expt%d%sAutoClusterTimesDetails.mat',fileparts(cfile),ex,xs);
    if exist(cfile)
        Clusters = LoadCluster(cfile,'rawxy');
    else
        Clusters = {};
    end
    if exist(afile)
        a = load(afile);
        if exist('ClusterDetails')  %have both
            for j = 1:length(a.ClusterDetails)
                if j > length(ClusterDetails) || isempty(ClusterDetails{j})
                    ClusterDetails{j} = a.ClusterDetails{j};
                end
            end
        else
            ClusterDetails = a.ClusterDetails;
        end
    end
    cfile = sprintf('%s/Expt%d%sClusterTimes.mat',fileparts(cfile),ex,xs);
end
if length(Clusters) >= probe
    Cp = Clusters{probe};
else
    Cp = {};
end
if length(Clusters) >= probe && isfield(Clusters{probe},'xy')
    Spikes.xy = Clusters{probe}.xy;
    Clusters{probe}.exptid = floor(Clusters{probe}.exptno);
elseif exist('ClusterDetails') && isfield(ClusterDetails{probe},'xy')
    Spikes.xy = Clusters{probe}.xy;
    if exist(cfile)
        load(cfile);
        Clusters{probe}.clst = ClusterDetails{probe}.clst;
        Clusters{probe}.xy = Spikes.xy;
    end
end

if isfield(Spikes,'xy')
    DATA.Spikes.cx = Spikes.xy(:,1);
    DATA.Spikes.cy = Spikes.xy(:,2);
end

if length(AllSpikes)
    AllSpikes{probes(1)}.xy = Spikes.xy;
    Spks = AllSpikes{probes(1)};
else
    Spks = Spikes;
end
if ~isfield(Spks,'values') %failed to load - 
    return;
end
if isfield(Spks,'xy') && size(Spks.xy,1) < size(Spks.values,1)
    fprintf('Size mismatch %d XY vals, %d Spike Vals\n',size(Spks.xy,1) ,size(Spks.values,1));
end
if isfield(C,'clst') && size(C.clst,1) == size(Spks.codes,1)
    Spks.codes(:,1) = C.clst -1;
    Spks.codes(:,2) = C.clst -1;
    if length(AllSpikes)
     AllSpikes{probes(1)}.codes(:,1) = C.clst-1;        
     AllSpikes{probes(1)}.codes(:,2) = C.clst-1;        
    end
end
    
t = Spks.times(1);
nt = 1;
tw = 10000;
if isempty(Expt)
    while t < max(Spks.times)
        Expt.Trials(nt).Start = t;
        Expt.Trials(nt).End = t + tw;
        if ~isfield(Expt.Trials,'Trial') || isempty(Expt.Trials(nt).Trials)
            Expt.Trials(nt).Trial = nt;
        end
        if ~isfield(Expt.Trials,'id') || isempty(Expt.Trials(nt).id)
            Expt.Trials(nt).id = nt;
        end
        t = t+tw+preperiod+postperiod;
        id = find(Spks.times > t);
        if length(id)
            t = Spks.times(id(1))-preperiod;
        end
        nt = nt+1;
    end
else
    if isfield(Expt.Header,'suffixes')
        sid = find(Expt.Header.suffixes ==ex);
        Expt.Header.BlockStartid(end+1) = 1+Expt.Trials(end).id;
        if ~isempty(sid)
            tid = find(ismember([Expt.Trials.id],Expt.Header.BlockStartid(sid):Expt.Header.BlockStartid(sid+1)));
            Expt.Trials = Expt.Trials(tid);            
        end
    end
    nt = length(Expt.Trials);
end
Expt.Header.trange = minmax(Spks.times);
Expt.Header.Name = spkfile;
Expt.Header.cellnumber = 0;
if ~isempty(parentfigure)
    DATA.toplevel = parentfigure;
end
DATA.plot.syncoverlay = 0;
if length(AllSpikes)
    DATA.AllSpikes = AllSpikes;
    DATA.plot.syncoverlay = 1;
else
DATA.AllData.Spikes = Spks;
end
DATA.Expts{1} = Expt;
DATA.spklist = 1:length(Spks.times);
DATA.AllData.pcs = [];
DATA.tag.clusterxy = 'Spike XY Plot';
DATA.tag.spikev = 'Spike Waveforms';
DATA.tag.top = DATA.tag.clusterxy;
DATA.state.uselfp = 0;
DATA.plot.clusterX = 1;
DATA.plot.clusterY = 2;
DATA.plot.clusterXrange = [0 40];
DATA.plot.clusterYrange = [0 5];
DATA.plot.clusterZrange = [0 5];
DATA.plot.clusterZ = 3;
DATA.plot.SpikeMaxV = max(abs(Spks.values(:)));
DATA.plot.SpikeVsep = 2;
DATA.plot.dvdt = 0;
DATA.plot.showwave = 0;
DATA.plot.prettyfigs = 1;
DATA.plot.synccluster = 0;
DATA.probe = probe;
DATA.currenttrial = 1;
DATA.plot.showsync = 0;
DATA.state.fixrange = 0;
DATA.state.preperiod = preperiod;
DATA.state.postperiod = postperiod;
DATA.syncsign = 6;
DATA.state.recut = 1;
DATA.cluster = [];
DATA.currentcluster = 1;
DATA.firsttrial = 1;
DATA.clusterArange = [7:11];
DATA.clusterBrange = [12:18];
DATA.clusterErange = [1:40];
DATA.plot.nodc = 1;
DATA.plot.voffsets(1) = 0; 
if length(probes) > 1
DATA.xprobes = probes(2:end);
np = (length(files)-1);
for j = 1:length(probes)
    range(j,1) = prctile(DATA.AllSpikes{probes(j)}.values(:),1);
    range(j,2) = prctile(DATA.AllSpikes{probes(j)}.values(:),99);
end
id = find(probes > probes(1));
id = sort(probes(id));
for j = 1:length(id);
    p = find(probes == id(j));
    q = find(probes == id(j)-1);
    DATA.plot.voffsets(p) = DATA.plot.voffsets(q)+range(p-1,2)-range(p,1);
end
id = sort(find(probes < probes(1)));
id = sort(probes(id));
for j = 1:length(id);
    p = find(probes == id(j)); %this probe
    q = find(probes == id(j)+1); %the probe above
    DATA.plot.voffsets(p) = DATA.plot.voffsets(q)+range(q,1)-range(q,2);
end
DATA.plot.voffsets = DATA.plot.voffsets .* 1.2;
end
    colors = mycolors;
    DATA.spkcolor{1} = [0.5 0.5 0.5];
    DATA.spkcolor(2:20) = colors(1:19);
    DATA.ptsize = 1;
    DATA.plot.xcorr = 0;
DATA.probelist = probes;
DATA.plot.autoscale = 1;
%DATA.AllData.Trialids = [Expt.Trials.Trial];
%DATA.Expts{1} = Expt;
DATA.explabels{1} = spkfile;
DATA.currentexpt = 1;

DATA.plot.showwave = 0;
DATA.vstep = 1;
DATA.densityplot = 0;
DATA.plot.setptsize = 0;
DATA.spikelist = [0:4];
DATA.sids = {};
for j = 1:length(excludetrials)
    DATA.Expts{1}.Trials(excludetrials(j)).Trial  =  abs(DATA.Expts{1}.Trials(excludetrials(j)).Trial)  .* -1; 
end


DATA = PC.SetDefaults(DATA);
DATA.mahaltype = distancemeasure;
DATA.options.usesavedcodes = 1;
DATA.nprobes = 24;
DATA.plot.density = 0;
if ~isfield(Spks,'Vrange')
    if isfield(Spks,'maxv')
        Spks.VRange = [-Spks.maxv Spks.maxv];
    else
        Spks.VRange = minmax(Spks.values(:));
    end
end
DATA.exptid = ex;
DATA.show.linecontextmenus=0;
DATA.plot.trighist = 1;

result.Clusters = Clusters;
result.probe = probe;
if quickspks
    [F, isnew] = GetFigure(DATA.tag.spikev,'parentfigure',parentfigure);
    if strcmp(DATA.Expts{1}.Header.DataType,'Spike2Swatches') % && selecttrial
        [tt, ts] = expt.BlockTimes(DATA.Expts{1});
        id = find(ts(:,1) < selecttrial & ts(:,2) > selecttrial);
        tt = tt .* 10000;
        sid = find(Spks.times > tt(id,1) & Spks.times < tt(id,2));
        Spks.times = Spks.times(sid);
        Spks.values = Spks.values(sid,:);
        Spks.codes = Spks.codes(sid,:);
    end
    PC.QuickSpikes(DATA,Spks,Cp);
    [F, isnew] = GetFigure(DATA.tag.clusterxy,'parentfigure',parentfigure);
    hold off;
    PC.PlotClusterXY(DATA,Cp);
else
    DATA = SpoolSpikes(DATA, varargin{:});
    result.toplevel = DATA.toplevel;
    result.svfig = DATA.svfig;
end

if isfield(DATA,'dprime')
    result.dprime = DATA.dprime;
end
if isfield(DATA,'dprimes')
    result.dprimes = DATA.dprimes;
end
if isfield(DATA,'dprimet')
    result.dprimet = DATA.dprimet;
end


function [probe, ex, xs] = Name2Probe(name)

cfile = strrep(name,'/Spikes','');
id = regexp(cfile,'p[0-9]*t[0-9]*');
if length(id)
probe = sscanf(cfile(id+1:end),'%d');
j = id+2;
while cfile(j) ~= 't'
    j = j+1;
end
ex = sscanf(cfile(j+1:end),'%d');
id = regexp(cfile,'p[0-9]*t[0-9]*a');

if length(id)
    xs = 'a';
else
    xs = [];
end
end


function Spikes = ReadSpikeFile(spkfile)
if ~exist(spkfile)
    fprintf('Cannot read %s\n',spkfile);
    Spikes = [];
    return;
end
X = load(spkfile);
if isfield(X,'Spikes')
    Spikes = X.Spikes;
else %old spike2 swatch file
    fprintf('Setting %s\n',spkfile);
    f = fieldnames(X);
    a = find(strncmp('Ch',f,2));
    Spikes = X.(f{a(1)});
    p = sscanf(Spikes.title,'Spike %d');
    clfile = strrep(spkfile,'/Spikes/','/');
    clfile = regexprep(clfile,'A.p[0-9]+t[0-9]+',sprintf('.p%dcl',p));
    clfile = regexprep(clfile,'.p[0-9]+t[0-9]+',sprintf('.p%dcl',p));
    if exist(clfile)
        X = load(clfile);
        clid  = X.clid;
    end
    probefile = regexprep(clfile,'.p[0-9]+cl','probes');
    X = load(probefile);    
    id = find([X.probes.probe] == p);
    [a,b] = fileparts(spkfile);
    fid = find(strcmp([b '.mat'],{X.probes(id).filename}));
    a = X.probes(id(fid)).firsti;
    Spikes.codes(:,1) = clid(a:a+Spikes.length-1);
    Spikes.codes(:,2) = clid(a:a+Spikes.length-1);
    Spikes.probe = p;
end
if isfield(Spikes,'Vrange') && ~isfield(Spikes,'VRange')
    Spikes.VRange = Spikes.Vrange;
end
Spikes.Header.loadname = spkfile;
if size(Spikes.values,2) > 100
    Spikes.values = Spikes.values';
end
    Spikes.times = Spikes.times .* 10000;
    Spikes.times = reshape(Spikes.times,length(Spikes.times),1);
    if size(Spikes.codes,2) == 1
        Spikes.codes(:,2) = Spikes.codes(:,1);
    end

