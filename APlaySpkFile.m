function [Expt, Expts, AllData, Raw] = APlaySpkFile(name, varargin)
% Lists Expts and Trials in .mat files from spike 2
% [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%
% Builds/retnes lists of trials/expts from a matlab file
% APlaySpkFile(name, 'relist')  Re-builds the index filt
% APlaySpkFile(name, 'setprobe', n)  sets which probes data is loaded
% initially
% APlaySpkFile(name, 'setprobe', -1) rebuilds probe list
% APlaySpkFile(name, 'usealltrials') includes bad fix trials too
% APlaySpkFile(name, 'expts') Just loads Expts file if availableo
%
% APlaySpkFile(name, 'rfs') rebulds .ufl file which record penetration
% info
%
% APlaySpkFile(..., 'makebnc') forces re-reading of a any .bnc file (online text
% made by manual expts). 

%Errors:  Can add txt lines post hoc in a file nameAdd.txt (i.e. jbeG044.mat -> jbeG044Add.txt)
%If it exists is read in and adds lines to the text stream.
%A file nameStimon  can be used to replace missing stimon/off DIO pulses
%A file nameAddEvents (EG jbeM016.14.AddEvents) can be used to add event codes
%  Format  time(seconds)  event number 
%
SpkDefs;
playspikes = 0;
defaults.fz = 96;
starttrial = 0;
onlinedata = 0;
idxfile = [];
Expts = [];
Expt = [];
AllData = [];
s2version = 0;
testframe = 0;
findtrial = 0;
spkch = 'Ch5';
stimch = 'Ch8';
framechname = 'Ch7';
mainsname = 'Ch24';
dstimch = 'Ch18';
ustimmarkname = '';
ustimmarkch = [];
mainsch = [];
stimlvl = [];
stimchange = [];
logfid = 0;
exptsonly = 0;
fixup = 1;
UstimV = 0;
setprobe = 1;
nerr = 0;
preperiod=2000;
postperiod = 2000;
%by default, do not set Stimulus.(et) etc from exvals.
%if this is neededed, set useexvals from Add.txt
useexvals = [0 0 0];
savedvdt = 0;
dvfile = [];
Raw = [];
timeoffset = 0;
ignoreSpikeO = 0;  %% for files with full voltage records elsewhere, ignore these
quickload = 0;
%method for matching up StimChan starts with text events. With many channel
%recordings, can get big delays between these, so the original method,
%finding the nearest stimchan event, is unreliable. New method (1) finds
%all events firts, and if the numbers are the same, matches in order.
%dfeault set to method 1 Jan 6 2012 by bgc
state.method = 1; 
state.showerrs = 2; %if 1, then popup windows are modal
saveexpts = 0;
state.needframes = 1;
state.tt = [];
state.nospikes = 0;
state.alltrials = 0;
state.profiling = 0;
state.resort = 0; %generate expt liat again. Else read from disk if there
state.buildbnc = 0;
state.autoinit = 0;
state.online = 0;
state.combineexpts = 0;
state.needbnc = 1;
state.bysuffix = 0;
state.keepspkvalues = 0;
verstr = '$Revision: 1.17 $';
state.version = sscanf(verstr(12:end),'%f');

errs = {};

mkidx =0;  %%need this up here so taht relist works
unsavedexpts = 0;

if isdir(name)
    [Expts, ~, Expt] = ReadExptDir(name, varargin{:});
    return; 
end

j = 1;
while j <= length(varargin)
    vg = varargin{j};
    if iscellstr(vg) %must be for someone else. Can't parse cell arrays here
    elseif strncmpi(vg,'alltrials',4)
        onlinedata = 1;
    elseif ischar(vg) & strncmpi(vg,'expts',4)
        exptsonly = 1;
    elseif strncmpi(vg,'Defaults',4)
        j = j+1;
        defaults = varargin{j};
        if isfield(defaults,'starttrial')
            starttrial = defaults.starttrial;
        end
        state = CopyFields(state,defaults,{'keepspkvalues'});
    elseif ischar(vg) & strncmpi(vg,'online',4)
        onlinedata = 1;
        state.online = 1;
        state.needbnc = 0;
    elseif ischar(vg) & strncmpi(vg,'needbnc',7)
        state.needbnc =1;

    elseif ischar(vg) & strncmpi(vg,'setprobe',4)
        j = j+1;
        setprobe = varargin{j};
    elseif strncmpi(vg,'findprobes',6)
       mkidx = 1;
       setprobe = -1; %force relisting of probes
    elseif strncmpi(vg,'modalerrs',5)
        state.showerrs = 1;
    elseif strncmpi(vg,'noerrs',5)
        state.showerrs = 0;
    elseif strncmpi(vg,'autoinit',5)
        state.autoinit = 1;
    elseif strncmpi(vg,'nospikes',5)
        state.nospikes = 1;
    elseif strncmpi(vg,'rfs',3)
       state.nospikes = 2;
    elseif strncmpi(vg,'bysuffix',7)
       ignoreSpikeO = 2;
    elseif strncmpi(vg,'method',6)
        j = j+1;
        state.method = varargin{j};
    elseif strncmpi(vg,'noframes',6)
        state.needframes = 0;

    elseif strncmpi(vg,'profile',6)
        state.profiling = 1;
    elseif sum(strncmpi(vg,{'quicksuffix' 'quickload'},9))
       ignoreSpikeO = 2;
       quickload = 1;
       if strcmp(vg,'quicksuffix')
           quickload = 2;
       end
    elseif strncmpi(vg,'quick',5)
        exptsonly = 1;
    elseif strncmpi(vg,'makebnc',4)
        state.buildbnc  = 1;
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    elseif strncmpi(vg,'resort',8)
        state.resort = 1;
    elseif strncmpi(vg,'saveexpts',6)
        saveexpts = 1;
    elseif strncmpi(vg,'sortexpts',6)
        idxfile = strrep(name,'.mat','idx.mat');
        load(idxfile);
        [Expts, Expt] = SortExpts(Expts, Expt.Trials, Expt.Header,1, Expt, state);
        return;
    elseif strncmpi(vg,'usealltrials',8)
        state.alltrials = 1;
    elseif strncmpi(vg,'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
    end
    
    j = j+1;
end


clusterdate = now;
thecluster = 1;
state.tt = TimeMark(state.tt,'Start');


if ischar(name)
    if ~exist(name,'file')
        fprintf('No file %s\n',name);
        return;
    end
    
    idxfile = strrep(name,'.mat','idx.mat');
    outname = strrep(name,'.mat','Expts.mat');
    if exptsonly && exist(outname);
        fprintf('Using %s\n',outname);
        [Expt, Expts] = LoadExptsFile(outname);
        if ~isfield(Expt,'Trials') & exist(idxfile)
            a = load(idxfile);
            if isfield(a,'Expt')
                Expt = a.Expt;
            end
        end
        return;
    end
    
    if strfind(name,'idx.mat')
        argon = {};
        load(name);
        if ~isfield(Expt,'state') 
            Expt.state.nospikes = 0;
        end
        Expt.Header.loadname = strrep(name,'idx.mat','.mat');
        state.tt = TimeMark(state.tt,'Loaded');
        Expt.state = CopyFields(Expt.state,state,{'showerrs'});
        [Expts, Expt, state] = SortExpts(Expts, Expt.Trials, Expt.Header, thecluster, Expt, state, argon{:});
        state.tt = TimeMark(state.tt,'Sorted');
        TimeMark(state.tt, state.profiling);
        return;
    end
    np = 0;
    nlfp = 0;
    nspkt = 0;
    probes = [];
    logname = strrep(name,'.mat', '.log');
%    fprintf('Log %s\n',logname);
   logfid = fopen(logname,'a');
   Oprobe = 0;
   AllSpks = {};
   if onlinedata
       %        oname = strrep(name,'online','online2');
       oname = strrep(name,'/Expt','A/Expt');
       oname = regexprep(oname,'(\.[0-9]*.mat)','A$1');
       idxfile = strrep(name,'.mat','idx.mat');
       if exist(idxfile,'file') && (~mkidx && state.keepspkvalues == 0)
           load(idxfile);
           probes = Expt.Probes;
       else
           if exist(oname,'file')
               af = load(oname);
               f = fields(af);
               if state.autoinit
                   AllData.Spike2ChA = af;
                   AllData.Spike2ChA.filename = oname;
                   AllData.Spike2ChA.source = 2;;
               end
               
           else
               f = {};
           end
       for j = 1:length(f)
           if ~isempty(regexp(f{j},'Ch[0-9]*'))
               ch = af.(f{j});
               if strncmpi(ch.title,'Spike',5)
                   np = np+1;
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                   if isempty(probe)
                       probes(np).probe = sscanf(vars{j},'Ch%d');
                   else
                       probes(np).probe = probe;
                   end
                   probes(np).var = f{j};
                   probes(np).traces = ch.traces;
                   probes(np).source = 2;
                   AllSpks{probes(np).probe}.values = ch.values;
                   AllSpks{probes(np).probe}.times = ch.times .* 10000;

                   if probe == setprobe(1)
                       Chspk = ch;
                       spkch = 'Chspk';
                   end
               elseif strncmpi(ch.title,'4Trode',5)
                   np = np+1;
                   probe = sscanf(ch.title,'4Trode%d');
                   if isempty(probe)
                       probes(np).probe = sscanf(vars{j},'Ch%d');
                   else
                       probes(np).probe = probe;
                   end
                   probes(np).var = f{j};
                   probes(np).traces = ch.traces;
                   probes(np).source = 2;
                   if probe == setprobe(1)
                       Chspk = ch;
                       spkch = 'Chspk';
                   end
               end
           end
       end
       end
   elseif ~isempty(regexp(name,'[0-9]\.[0-9]*\.mat'))
       state.bysuffix = 1;
       if state.resort || mkidx || state.online == 0
           state.nospikes = 1;
       end
   end
%    fprintf('Reading %s\n',name);
   mkmatver = 0;


   if ignoreSpikeO == NaN % don't need this any more, for online at least
       oname = regexprep(name,'.([0-9]*.mat)','A.$1');
       if exist(oname,'file')
           load(oname);
       end
       avars = who('Ch[0-9]*');
       for j = 1:length(avars)
           eval([avars{j} 'A = ' avars{j} ';']);
           clear(avars{j});
       end

   end

   ts = now;
   load(name);
   Expt.loaddur(1) = mytoc(ts);
    
    if exist('SMRFiles','var') %% Concatentae existing files
        toff = 0;
        for j = 1:length(SMRFiles.Names)
            [a,b,c] = APlaySpkFile(SMRFiles.Names{j},varargin{:},'timeoffset',toff);
            if j == 1
                Expt = a;
                Expts = b;
                AllData = c;
            else
                f = fields(Expt.Trials);
                for k = 1:length(f)
                    if isfield(a.Trials,f{k})
                        if diff(size(Expt.Trials.(f{k}))) > 1
                            dim = 2;
                        else
                            dim = 1;
                        end
                    Expt.Trials.(f{k}) = cat(dim,Expt.Trials.(f{k}), a.Trials.(f{k}));
                    end
                end
                Expt.Probes = [Expt.Probes a.Probes];
                Expt.Spkid = cat(1,Expt.Spkid,a.Spkid);
                Expts = [Expts b];
                f = fields(AllData.Spikes);
                for k = 1:length(f)
                AllData.Spikes.(f{k}) = cat(1,AllData.Spikes.(f{k}),c.Spikes.(f{k}));
                end
            end
            toff = max(Expt.Trials.End)+1;

        end
        return;
    end
    if state.nospikes == 2  %make ufl file
        MkUfl(name,Ch30,'overwrite');
        return;
    end
    if exist('Ch31','var') & isfield(Ch31,'comment')
        mkmatver = sscanf(Ch31.comment,'MkMat V%n');
    end
    if exist('Ch30','var') & strncmp(Ch30.comment,'GridData',8) 
        state.nospikes = 1; %Don't try to load up all spike files if its Utah Array
    end
    
    
    Spkerrs = [];;
    vars = who('Ch*');
    xv = who('Failed*');
    if ~isempty(xv)
        cprintf('red','File %s has export failues\n',name);
    end
    vars = {vars{:} xv{:}};
    for j = 1:length(vars)
        if state.autoinit
            AllData.Spike2Ch.(vars{j}) = eval(vars{j});
            AllData.Spike2Ch.filename = name;
            AllData.Spike2Ch.source = 1;
        end
        if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*')) || strncmp(vars{j},'Failed',5)
            eval(['ch = ' vars{j} ';']);
            chn = sscanf(vars{j},'Ch%d');
            if chn > 400  %a memory/extra offline channel
            elseif strncmpi(ch.title,'SpikeO',6) && ignoreSpikeO;
            elseif strncmpi(ch.title,'Spike',5) && state.nospikes == 0
                np = np+1;  
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
                    probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                probes(np).traces = ch.traces;
                probes(np).source = 1;
                AllSpks{probes(np).probe}.values = ch.values;
                AllSpks{probes(np).probe}.times = ch.times .* 10000;


            elseif strncmpi(ch.title,'4Trode',5)
                np = np+1;  
                probe = sscanf(ch.title,'4Trode%d');
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                if probe == setprobe(1)
                    Chspk = ch;
                    spkch = 'Chspk';
                end
                probes(np).traces = ch.traces;
                probes(np).source =1;

            elseif strncmpi(ch.title,'uStimMk',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            elseif strncmpi(ch.title,'uStim',5)
                np = np+1;  
                probe = 100;
                UstimV = ch;
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
                    probes(np).probe = probe;
                end
                probes(np).source = 0;
                probes(np).traces = ch.traces;
                probes(np).var = vars{j};
                
            elseif strncmpi(ch.title,'StimOn',6)
                stimch = vars{j};
                stimlvl = ch;
            elseif strncmpi(ch.title,'StimChan',8) % stim change detector
                dstimch = vars{j};
                stimchange = ch;
            elseif strncmpi(ch.title,'VTR',3)
                framechname = vars{j};
                framech = ch;
   %             fprintf('Frames in %s\n',vars{j});
            elseif strncmpi(ch.title,'Mains',5)
                mainsname = vars{j};
                mainsch = ch;
            elseif strncmpi(ch.title,'DigMark',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            end
        end
    end
    
    
    if np > 0
       if np == 1
           probes(1).probe = 1;
       end
        if exist('Chspk','var') %found selected probe
            spkch = 'Chspk';
        elseif setprobe(1) > 0 
            id = find([probes.probe] == setprobe(1));
            if ~isempty(id)
             spkch = probes(id(1)).var;
            else
             spkch = probes(1).var;
            end
        else
        spkch = probes(1).var;
        end
    else
        idxfile = strrep(name,'.mat','probes.mat');
        if exist(idxfile,'file') & setprobe(1) >= 0
            load(idxfile);
            id = find([probes.probe] == setprobe(1));
            [dp, pref] = fileparts(name);
            spkdir = [dp '/Spikes'];
 % if IgnoreSpikeO is 2 this means don't load the spike waveforms, just the cluster
 % cluster xy values from the ClusterDetails File
            if ~exist(spkdir,'dir') || ignoreSpikeO == 2
                spkdir = dp;
            end

            
            if length(id) > 1 
                ch = sscanf(probes(id(1)).var,'Ch%d');
            else
                ch = sscanf(probes(1).var,'Ch%d');
            end

            if isfield(probes,'filename')
                pname = [spkdir '/' probes(id(1)).filename];
            elseif setprobe > 16
               pname = strrep(name,'.mat',sprintf('A.p%d.mat',ch));
            else
               pname = strrep(name,'.mat',sprintf('.p%d.mat',ch));
            end
%
%  Need to change this to load multiple file if necessary

            dvfile = [spkdir '/' strrep(probes(id(1)).filename,'t0.','dvdt.')];
            if length(id) > 1 & length(setprobe) == 1
                [a,sid] = sort([probes(id).first]);
                sid = id(sid);
                Spk.times = [];
                Spk.values = [];
                Spk.codes = [];

                for j = 1:length(sid)
                    filename = [dp '/Spikes/' probes(sid(j)).filename];
                    if exist(filename,'file')
                        a = load(filename);
                        chname = probes(sid(j)).var;
                        if isempty(Spk.times)
                            Spk = a.(chname);
                        else
                            Spk.times = [Spk.times; a.(chname).times];
                            Spk.codes = [Spk.codes; a.(chname).codes];
                            Spk.values = [Spk.values; a.(chname).values];
                            if size(a.(chname).times,1) > size(a.(chname).values,1)
                                fprintf('Some Missing Spike values in %s\n',filename);
                                if logfid
                                    fprintf(logfid,'Some Missing Spike values in %s\n',filename);
                                end
                            end
                        end
                    else
                        fprintf('No file %s\n',filename);
                    end
                end
                spkch = 'Spk';
            else
                a = load(pname);
                chid = find(strncmp('Ch',fieldnames(a),2));
                chnames = fieldnames(a);
                Chspk = a.(chnames{chid(1)});
                spkch = 'Chspk';
            end
                clear('a');

        elseif state.nospikes == 0
            if setprobe < 0 %forces rebuilding of probe list
                setprobe = -setprobe;
            end
            [a,b] = splitpath(name);
            if ignoreSpikeO == 2
                d = dir([b '/Expt*Times.mat']);
            else
                d = dir(b);
            end
                
            spkdir = [b '/Spikes'];
            if exist(spkdir,'dir') && ignoreSpikeO ~= 2
                d = dir(spkdir);
            else
                spkdir = b;
            end
            Chspk = [];
            for j = 1:length(d)
                if regexp(d(j).name,'.p[0-9][0-9,t]*.mat') %% d(j).bytes < 525000000 & np < 2
                    fprintf('loading %s\n',d(j).name)
                    a =  load([spkdir '/' d(j).name]);
                    chnames = fieldnames(a);
                    chnames = chnames(find(strncmp('Ch',chnames,2)));
                    if ~isempty(chnames) %is data in file
                        np = np+1;
                        if strncmp(a.(chnames{1}).title,'4Trode',6)
                            probes(np).probe = sscanf(a.(chnames{1}).title,'4Trode%d');
                        else
                            probes(np).probe = sscanf(a.(chnames{1}).title,'Spike %d');
                        end
                        probes(np).var = chnames{1};
                        probes(np).filename = d(j).name;
                        probes(np).pathname = spkdir;
                        probes(np).first = a.(chnames{1}).times(1);
                        probes(np).last = a.(chnames{1}).times(end);
                        probes(np).nspk = length(a.(chnames{1}).times);
                        probes(np).traces = a.(chnames{1}).traces;
                        if strfind(d(j).name,'A.p')
                        probes(np).source = 2;
                        else
                        probes(np).source = 1;
                        end
                        if probes(np).probe == setprobe
                            Chspk = a.(chnames{1});
                        end
                        if length(a.(chnames{1}).times) > size(a.(chnames{1}).values,1)
                            nerr=nerr+1;
                            errs{nerr} = [d(j).name 'Missing some values'];
                            if state.showerrs == 1
                                msgbox(errs{nerr},'APlaySpkFile Error!!','modal');
                            end
                            fprintf('%s\n',errs{nerr});
                        end
                    end
                elseif regexp(d(j).name,'.p[0-9][0-9]*.mat') %% d(j).bytes < 525000000 & np < 2
                    fprintf('loading %s\n',d(j).name)
                    a =  load([spkdir '/' d(j).name]);
                    chnames = fieldnames(a);
                    chnames = chnames(find(strncmp('Ch',chnames,2)));
                    if ~isempty(chnames) %is data in file
                        np = np+1;
                        probes(np).probe = sscanf(a.(chnames{1}).title,'Spike %d');
                        probes(np).var = chnames{1};
                        probes(np).filename = d(j).name;
                        probes(np).first = a.(chnames{1}).times(1);
                        probes(np).last = a.(chnames{1}).times(end);
                        if probes(np).probe == setprobe
                            Chspk = a.(chnames{1});
                        end
                        if length(a.(chnames{1}).times) > size(a.(chnames{1}).values,1)
                            nerr=nerr+1;
                            errs{nerr} = [d(j).name 'Missing some values'];
                            msgbox(errs{nerr},'APlaySpkFile Error!!','modal');
                            fprintf('%s\n',errs{nerr});
                        end
                    end
                elseif regexp(d(j).name,'NewClusterTimes.mat') %% Ignore these
                elseif regexp(d(j).name,'ClusterTimes.mat') %% d(j).bytes < 525000000 & np < 2
                    suffix = 0;
                    csuffix = 0;
                    id = regexp(name,'\.[0-9,a]*\.mat');
                    if length(id)
                        suffix = sscanf(name(id(1)+1:end),'%d');
                        id = regexp(d(j).name,'Expt.[0-9,a]*');
                        if length(id)
                            csuffix= sscanf(d(j).name(id(1)+4:end),'%d');
                        end
                    end
%shouldn't need to load autoclustertimes except special circumstances
%since its put into ClusterTimes. When would you combine without making
%ClusterTimes? 
% - when the autoclusters are good!  Forget this can happen. 
%Only load autoclustertimes if clusertimes doesn't exist
                    if ~isempty(strfind((d(j).name),'AutoCluster'))
                        s = strrep(d(j).name,'AutoCluster','Cluster');
                        if sum(strcmp(s,{d.name}))
                            go = 0;
                        else
                            go = 1;
                        end
                    else
                        go = 1;
                    end
                    if quickload == 2 %temp test for dynamically loading everything
                       go = 0;
                    end
                    if suffix == csuffix && go
                    fprintf('loading %s\n',d(j).name);
                    clusterdate = d(j).datenum;
                    cname = [spkdir '/' d(j).name];
                    te = now;
                    a =  load(cname);
                    a.filename = cname;
                    if quickload == 0
                        dname = strrep(cname,'ClusterTimes','ClusterTimesDetails');
                        if exist(dname)
                            b = load(dname);
                            gotempty = 0;
                            for k = 1:length(b.ClusterDetails)
                                if isfield(b.ClusterDetails{k},'t')
                                    a.Clusters{k}.times = b.ClusterDetails{k}.t;
                                    a.Clusters{k}.xy = b.ClusterDetails{k}.xy;
                                    a.Clusters{k}.clst = b.ClusterDetails{k}.clst;
                                    %                            a.Clusters{k}.Evec = b.ClusterDetails{k}.Evec;
                                    for c = 1:size(b.ClusterDetails{k}.next);
                                        if isfield(b.ClusterDetails{k}.next{c},'xy')
                                            a.Clusters{k}.next{c}.xy = b.ClusterDetails{k}.next{c}.xy;
                                        end
                                    end
                                else
                                    gotempty = gotempty+1;
                                end
                            end
                        else
                            gotempty = 1;
                        end
                        if gotempty || length(b.ClusterDetails) < length(a.Clusters)
                            dname = strrep(cname,'ClusterTimes','AutoClusterTimesDetails');
                            if exist(dname)
                                b = load(dname);
                                for k = 1:length(a.Clusters)
                                    if isfield(b.ClusterDetails{k},'t') && ~isfield(a.Clusters{k},'xy')
                                        a.Clusters{k}.times = b.ClusterDetails{k}.t;
                                        a.Clusters{k}.xy = b.ClusterDetails{k}.xy;
                                        a.Clusters{k}.clst = b.ClusterDetails{k}.clst;
                                        %                            a.Clusters{k}.Evec = b.ClusterDetails{k}.Evec;
                                    end
                                end
                            end
                        end
                    else
                    end
                    nspkt = nspkt+1;
                    SpkTimes(nspkt).loadtime = (now-te)*60 *60*24;
                    id = strfind(d(j).name,'Expt');
                    if id
                        SpkTimes(nspkt).expno = sscanf(d(j).name(id(1)+4:end),'%d');
                    else
                        SpkTimes(nspkt).expno = 0;
                    end
                    SpkTimes(nspkt).section = 0;
                    id = regexp(name,'\.[0-9]*a\.mat');
                    if length(id)
                        SpkTimes(nspkt).section = 1;
                    end
                    
                    if strfind((d(j).name),'AutoCluster')
                        AutoClusters = a.Clusters;
                        SpkTimes(nspkt).auto = 1;
                    else
                        SpkTimes(nspkt).auto = 0;
                    end
                    SpkTimes(nspkt).datenum = clusterdate;
                    SpkTimes(nspkt).filename = d(j).name;
                    SpkTimes(nspkt).loadname = a.filename;
                    mint = 10e12;
                    maxt = 0;
                    for k = 1:length(a.Clusters)
                        if ~isempty(a.Clusters{k})
                            if a.Clusters{k}.auto == 0 && SpkTimes(nspkt).auto == 0
                                a.Clusters{k}.auto = 0;
                            else
                                a.Clusters{k}.auto = 1;
                            end
                            mint = min([mint min(a.Clusters{k}.times)]);
                            maxt = max([mint max(a.Clusters{k}.times)]);
                        end
                        if quickload == 0 && ~isfield(a.Clusters{k},'xy') && isfield(SpkTimes(1),'Clusters') &&isfield(SpkTimes(1).Clusters{k},'xy')
                            a.Clusters{k}.xy = SpkTimes(1).Clusters{k}.xy;
                            a.Clusters{k}.times = SpkTimes(1).Clusters{k}.times;
                            a.Clusters{k}.clst = SpkTimes(1).Clusters{k}.clst;
                        end
                        if diff(size(a.Clusters{k}.times)) < 0 %make sure these are column vectors
                            cprintf('red','%s P%d has t as row vector\n',cname,k)
                            a.Clusters{k}.times = a.Clusters{k}.times';
                        end
                    end
                    SpkTimes(nspkt).Clusters = a.Clusters;
                    if isfield(a,'FullVData')
                        SpkTimes(nspkt).FullVData = a.FullVData;
                    elseif isfield(b,'FullVData')
                        SpkTimes(nspkt).FullVData = b.FullVData;
                    end
                    SpkTimes(nspkt).trange = [mint maxt];
                    end
                    
                elseif regexp(d(j).name,'ClusterTimesDetails.mat') %% d(j).bytes < 525000000 & np < 2
                end
            end
            if isempty(Chspk) & exist('chnames','var') && ~isempty(chnames)
                Chspk = a.(chnames{1});
            end
            if ~isempty(probes)
                [a, id] = sort([probes.probe]);
                probes = probes(id);
                ip = unique([probes.probe]);
%sort by time of first spike keep track of total # spikes
                for j = ip;
                    id = find([probes.probe] == j);
                    [t, sid] = sort([probes(id).first]);
                    probes(id) = probes(id(sid));
                    nspk = 1;
                    for k = id;
                        probes(k).firsti = nspk;
                        nspk = nspk + probes(k).nspk;
                    end
                end
                save(idxfile,'probes');
            end
            clear a;
            spkch = 'Chspk';
        else %state.nospikes > 0
            [a,b] = splitpath(name);
            d = dir([b '/*ClusterTimes.mat']);
        end

    end
    if ~isempty(probes)
        [a, id] = sort([probes.probe]);
        probes = probes(id);
    end
    vnames = {'Ch30' 'Ch31' stimch};
    vlabels = {'Text' 'Events' 'Stim ON/OFF'};
    if state.needframes
        vlabels = [vlabels 'Frames'];
        vnames = [vnames framechname];
    end
    if state.nospikes == 0
        vlabels = [vlabels 'Spikes'];
        vnames = [vnames spkch];
    end
    for j = 1:length(vnames)
        missing(j) = ~exist(vnames{j},'var');
    end
    if sum(missing)
        msgbox(sprintf('%s Missing %s',name,vlabels{find(missing)}),'APlaySpkFile Error!!','modal');
        fprintf('%s Missing %s\n',name,vlabels{find(missing)});
        if logfid
            fprintf(logfid, '%s Missing %s\n',name,vlabels{find(missing)});
        
        fclose(logfid);
        end
        return;
    end
    Text = Ch30;
    if state.nospikes
        Spks.loaded = 0;
        Spks.times = [];
        Spks.codes = [];
        Spks.probe = 0;
    elseif ignoreSpikeO == NaN && length(probes) > 1   && Oprobe > 1
        for j = 1:length(probes)
            spkch = probes(j).var;
            p = probes(j).probe;
            AllSpikes{p} = eval(['CleanSpikes(' spkch ');']);
            AllSpikes{p}.times = AllSpikes{p}.times + timeoffset./10000;
        end
        Spks = AllSpikes;
    elseif ~eval(['isempty(' spkch ')']);
        if savedvdt & dvfile
            Spks = eval(['CleanSpikes(' spkch ',''dvfile'',dvfile);']);
        else
            Spks = eval(['CleanSpikes(' spkch ');']);
        end
        Spks.times = Spks.times+timeoffset./10000;
        clear(spkch);
    else
        if nspkt
            expts = unique([SpkTimes.expno]);
            expts = expts(expts > 0);
            np = 0;
            for j = expts
                id = find([SpkTimes.expno] == j & [SpkTimes.auto] == 0);
                aid = find([SpkTimes.expno] == j & [SpkTimes.auto] == 1);
                if length(aid) & length(id)

                    Clusters = SpkTimes(id(1)).Clusters;
                    for c = 1:length(SpkTimes(aid(1)).Clusters);
                        if c > length(Clusters) || isempty(Clusters{c})
                            Clusters{c} = SpkTimes(aid(1)).Clusters{c};
                        end
                        if ~isfield(Clusters{c},'excludetrialids')
                            Clusters{c}.excludetrialids{1} = [];
                        else
                            c = c;
                        end
                    end
                    if length(aid) == 2 && length(id) ==2 %Expt was split up
                        bClusters = SpkTimes(id(2)).Clusters;
                        for c = 1:length(SpkTimes(aid(2)).Clusters);
                            if c > length(bClusters) || isempty(bClusters{c})
                                Clusters{c}.times  = cat(2,Clusters{c}.times, SpkTimes(aid(2)).Clusters{c}.times);
                                Clusters{c}.clst  = cat(1,Clusters{c}.clst, SpkTimes(aid(2)).Clusters{c}.clst);
                                Clusters{c}.xy  = cat(1,Clusters{c}.xy, SpkTimes(aid(2)).Clusters{c}.xy);
                            else
                                Clusters{c}.times  = cat(2,Clusters{c}.times, bClusters{c}.times);
                                Clusters{c}.clst  = cat(1,Clusters{c}.clst, bClusters{c}.clst);
                                Clusters{c}.xy  = cat(1,Clusters{c}.xy, bClusters{c}.xy);
                                if isfield(bClusters{c},'excludetrialids')
                                    Clusters{c}.excludetrialids{1} = cat(2,Clusters{c}.excludetrialids{1}, bClusters{c}.excludetrialids);
                                end
                                    
                            end
                        end
                    end
                elseif length(aid)
                    Clusters = SpkTimes(aid).Clusters;
                else
                    Clusters = SpkTimes(id).Clusters;
                end
                np = max([np length(Clusters)]);
                EClusters{j} = Clusters;
                suffixlist(j) = j;
            end
            for j = 1:np
                Spks{j}.times = [];
                Spks{j}.codes = [];
                Spks{j}.cx = [];
                Spks{j}.cy = [];
                Spks{j}.errs = {};
                Spks{j}.savetime = [];
            end
            for k = 1:length(EClusters)
                Clusters = EClusters{k};
            for j = 1:length(Clusters)
                if Clusters{j}.shape == 0
                    Clusters{j}.sign = 1;
                    Clusters{j}.crit = 0;
                end
                if ~isfield(Clusters{j},'sign') || Clusters{j}.sign == 0
                    Clusters{j}.sign = 1;
                end
                if isfield(Clusters{j},'savetime')
                    Spks{j}.savetime = Clusters{j}.savetime(end);
                end
                if isfield(Clusters{j},'times')
                    probes(j).probe = j;
                    if isfield(Clusters{j},'errs')
                        Spks{j}.errs = Clusters{j}.errs;
                        for k = 1:length(Clusters{j}.errs)
                            Spkerrs(end+1).s = sprintf('%s',deblank(Clusters{j}.errs{k}));
                            Spkerrs(end).p = Clusters{j}.probe(1);
                            Spkerrs(end).exptno = Clusters{j}.exptno;
                        end
                    else
                        Spks{j}.errs = {};
                    end
                    if isempty(Clusters{j}.times)
                        Spks{j}.times = [];
                    elseif diff(size(Clusters{j}.times) > 1)
                        Spks{j}.times = cat(1,Spks{j}.times,Clusters{j}.times');
                    else
                        Spks{j}.times = cat(1,Spks{j}.times,Clusters{j}.times);
                    end
                    if quickload %just loaded times of classified spikes
                        Spks{j}.codes(:,1) = ones(size(Clusters{j}.times));
                        for c = 1:length(Clusters{j}.next)
                            if isfield(Clusters{j}.next{c},'times')
                            x = [Spks{j}.codes(:,1)' ones(size(Clusters{j}.next{c}.times)) * (1+c)];
                            Spks{j}.codes = x';
                            Spks{j}.times = cat(1,Spks{j}.times,Clusters{j}.next{c}.times');
                            end
                        end
                    else
                        if isfield(Clusters{j},'clst')
                            Spks{j}.codes(:,1) = Clusters{j}.clst-1;
                        else
                            Spks{j}.codes(:,1) = Clusters{j}.xy(:,1).*Clusters{j}.sign > Clusters{j}.crit.*Clusters{j}.sign;
                        end
                        Spks{j}.cx = Clusters{j}.xy(:,1);
                        Spks{j}.cy = Clusters{j}.xy(:,2);
                    end
%                    Spks{j}.codes(:,2) = Spks{j}.codes(:,1);
%mahal 1 9s 2D, mahal 4 in 1D, mahal 3 id ND, 0 if not used
                    if isfield(Clusters{j},'fitdprime')
                        Spks{j}.dips = [Clusters{j}.fitdprime(1) Clusters{j}.mahal(4) Clusters{j}.mahal(1)];
                        if Clusters{j}.space(1) == 6 && Clusters{j}.mahal(3) > Clusters{j}.mahal(1)
                            Spks{j}.dips(3) = Clusters{j}.mahal(3);
                        end
                    elseif isfield(Clusters{j},'mahal')
                        Spks{j}.dips = [NaN Clusters{j}.mahal(4) Clusters{j}.mahal(1)];
                    elseif isfield(Clusters{j},'dipsize')
                        Spks{j}.dips = [Clusters{j}.hdip(1) Clusters{j}.dipsize(1) abs(Clusters{j}.dprime(1))];
                    else
                        Spks{j}.dips = [NaN NaN NaN];
                    end
                    if isfield(Clusters{j},'mydip')
                        Spks{j}.dips(4) = Clusters{j}.mydip(3);
                    else
                        fprintf('Missing mydip in %d\n',j);
                    end
                    if isfield(Clusters{j},'dropi')
                        Spks{j}.dropi = Clusters{j}.dropi(3);
                    end
                    Spks{j}.next = [];
                    if isfield(Clusters{j},'next') && iscell(Clusters{j}.next)
                        for k = 1:length(Clusters{j}.next)
                            if isfield(Clusters{j}.next{k},'fitdprime')
                                Spks{j}.next{k}.dips(1) = Clusters{j}.next{k}.fitdprime(1);
                            end
                            if isfield(Clusters{j}.next{k},'mahal')
                                Spks{j}.next{k}.dips(2) = Clusters{j}.next{k}.mahal(4);
                                Spks{j}.next{k}.dips(3) = Clusters{j}.next{k}.mahal(1);
                            end
                            if isfield(Clusters{j}.next{k},'dropi')
                                Spks{j}.next{k}.dropi = Clusters{j}.next{k}.dropi(3);
                            end
                            if isfield(Clusters{j}.next{k},'shape')
                                Spks{j}.next{k}.shape = Clusters{j}.next{k}.shape;
                                Spks{j}.next{k} = CopyFields(Spks{j}.next{k},Clusters{j}.next{k},{'xy' 'xyr' 'angle' 'shape'});
                            end
                        end
                    end
                    Spks{j}.mahal = [Clusters{j}.mahal(1) Clusters{j}.mahal(2)];
                    Spks{j}.suffix = GetExptNumber(SpkTimes.filename); 
                    Spks{j}.crit = Clusters{j}.crit;
                    Spks{j}.clustersign = Clusters{j}.sign;
                    if isfield(Clusters{j},'excludetrialids')
                        Spks{j}.excludetrialids{1} = Clusters{j}.excludetrialids;
                    else
                        Spks{j}.excludetrialids{1} = [];
                    end
                    if isfield(Clusters{j},'missingtrials')
                        Spks{j}.excludetrialids{1} = union(Spks{j}.excludetrialids{1},Clusters{j}.missingtrials);
                    end
                    if length(Spks{j}.excludetrialids{1}) > 0
                        a = max(Spks{j}.excludetrialids{1});
                    end
                end
            end
            end
        else
            Spks.times = 0;
            Spks.codes = [0 0 0 0];
            Spks.values = 0;
        end
    end
    Events = Ch31;
    idxfile = strrep(name,'.mat','idx.mat');
    if idxfile & ~exist(idxfile,'file')
        mkidx = 1;
    elseif mkidx < 0
        mkidx = 0;
    elseif idxfile % idx file must exist - check its date
            d = dir(idxfile);
            if isempty(d)
                fprintf('Dir is Missing exisitng index file %s!!!!!!!\n',name,idxfile);
                mkidx = 1;
            else
                dd = datenum(d.date);
                d = dir(name);
                if datenum(d.date) > dd % Expt file newer - rebuild
                    fprintf('Rebuilding index: %s newer than %s\n',name,idxfile);
                    mkidx = 1;
                end
            end
    end
  
end


if isstruct(name)
    if isfield(name,'title') & isfield(name,'values') &... 
            size(name.values,2) == 46
        Spks = name;
        name = 'Unnamed'
    end
end



stimlvl = AddStimsFromFile(strrep(name,'.mat','.Stimon'),stimlvl);
Events = AddEventsFromFile(strrep(name,'.mat','.AddEvents'),Events);
forcefix = 0;
argon = {};
j = 1;
while j <= nargin-1
    vg = varargin{j};
    if isstruct(vg) 
        if isfield(vg,'text')
            Text = vg;
        elseif isfield(vg,'codes');
            Events = vg;
        end
    elseif strncmpi(vg,'alltrials',8)
        argon = {argon{:} varargin{j}};
    elseif strncmpi(vg,'Defaults',4)
        j = j+1;
        defaults = varargin{j};
        if isfield(defaults,'starttrial')
            starttrial = defaults.starttrial;
        end
        state = CopyFields(state,defaults,{'keepspkvalues'});
    elseif strncmpi(vg,'fixlfp',5)
        if strncmpi(vg,'fixlfpforce',8)
            forcefix = 1;
        end
        lfpfile = strrep(name,'.mat','A.lfp.mat');
        fixfile = strrep(name,'.mat','.lfp.mat');
% with 8 channels, everything is in one file - there is no 'A.lfp.mat' 
        if ~exist(lfpfile,'file') && exist(fixfile,'file')
            lfpfile = fixfile;
            forcefix = 1;
        end
        if exist(lfpfile,'file') && (~exist(fixfile,'file') || forcefix) ...
                && exist(mainsname,'var')
            load(lfpfile);
            if isfield(LFP.Header,'MainsNoise')
                fprintf('%s Already fixed\n',lfpfile);
                if logfid
                fprintf(logfid,'%s Already fixed\n',lfpfile);
                end
            else
            [LFP, avgs, NoiseAmp] = FixLFPMains(LFP,mainsch.times .* 10000);
            LFP.Header.amps = LFPGains(LFP);
            a = LFP.Header.amps ./ max(LFP.Header.amps);
            nch = sum(LFP.Header.chanlist > 0);
            if std(a(find(a > 0.1))) > 0.2 && nch <= 8  %% 8 channel probe tends to have mixed LFP gain
                LFP.Header.needscale = 1;
            else
                LFP.Header.needscale = 0;
            end
            save(fixfile,'LFP','NoiseAmp');
           if logfid > 0
               fprintf(logfid, '%s Fixed LFP in %s\n',datestr(now),fixfile);
           end
            end
        else
            fprintf('No LFP file %s\n',lfpfile);
           if logfid > 0
               fprintf(logfid, 'No LFP file %s\n',lfpfile);
           end
            
        end
        if strncmpi(vg,'fixlfponly',8)
            return;
        end
            
    elseif strncmpi(vg,'findtrial',5)
        j = j+1;
        argon = {argon{:} varargin{j-1} varargin{j}};
        findtrial = varargin{j};
    elseif strncmpi(vg,'mkufl',4)
        MkUfl(name, Ch30);
    elseif strncmpi(vg,'name',4)
        j = j+1;
        name = varargin{j};
   elseif strncmpi(vg,'noidx',5)
       mkidx = 0;
    elseif strncmpi(vg,'online',4)
        onlinedata = 1;
    elseif strncmpi(vg,'cluster',2)
        j = j+1;
        thecluster = varargin{j};
    elseif strncmpi(vg,'play',4)
        playspikes = 1;
    elseif strncmpi(vg,'rfs',3)
        MkUfl(name, Text,'overwrite');
    elseif strncmpi(vg,'relistonly',8)
        mkidx = 2;
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    end
    j = j+1;
end



 if isfield(defaults,'toplevel')
    state.toplevel = defaults.toplevel;
 end
Header.Name = BuildName(name);
Header.loadname = name; 

Events.times = Events.times * 10000;
%spk times need to be ints for trigsdf.
if ~isempty(Spks)
    if iscell(Spks)
        for j = 1:length(Spks)
            if isfield(Spks{j},'times')
                Spks{j}.times = round(Spks{j}.times * 10000);
            end
        end
    else
        Spks.times = round(Spks.times * 10000);
    end
end
if isstruct(UstimV) && ~isempty(UstimV)
UstimV.times = round(UstimV.times * 10000);
end
Text.times = Text.times * 10000;

if strncmp(Text.comment,'GridData',8) || strncmp(Text.comment,'uProbe',8)
    Expt.DataType = Text.comment;
else
    Expt.DataType = 'Spike2';
end
Expt.progname = 'APlaySpkFile';
Expt.exptno = GetExptNumber(name);
Expt.Header.DataType = Expt.DataType;
Expt = ShowErrors(Expt, Spkerrs); %Add in errors from ClusterFiles

if ~mkidx
    load(idxfile);
    if ~isfield(Expt,'errs')
        Expt.errs = {};
    end
    if ~isfield(Expt,'DataType')
        Expt.DataType = 'Spike2';
    end
    if ~isfield(Expt,'progname')
       Expt.progname = 'APlaySpkFile';
       Expt.exptno = GetExptNumber(name);
    end
    Expt.newerrs = 0;
% make sure name is is correct windows/unix form
    Expt.Header.Name = BuildName(Expt.Header.Name);
    Expt.Header.loadname = Header.loadname;
    Expt.setprobe = setprobe;
    Expt.Header.DataType = Expt.DataType;
    aText.text = {};
    aText.times = [];
    aText.codes = [];
    [aText, Ch30, na] = AddText(regexprep(name,'\.[0-9]*.mat','Add.txt'), aText, Ch30);
    [aText, Ch30, nb] = AddText(strrep(name,'.mat','Add.txt'), aText, Ch30);
%only re-read peninfo if lines are added    
    if isfield(Expt,'Comments') && isfield(Expt.Comments,'Peninfo') && (na+nb) > 0
        if isempty(aText) || isempty(aText.text);
            txt = GetPenInfo(Ch30,'name',name);
        else
            txt = GetPenInfo(aText,'name',name);
        end
        txt = strrep(txt,'Contact CNT','ContactCNT');
        Expt.Comments.Peninfo.trode = txt;
        id = strfind(txt,'Contact');
        if length(id)
            x = id(1);
            id = strfind(txt(x:end),' ');
            x = sscanf(txt(id(1)+x:end),'%d');
            Expt.Comments.Peninfo.probesep = x;
        end
        [rfstr, rfdat] = MkUfl(name, Text);
        Expt.Header.rfstr = rfstr;
        Expt.Header.rf = rfdat;
    end
    if iscell(Spks)  && isfield(Spks{1},'values')
        AllData.AllSpikes = Spks;
        Spks = AllData.AllSpikes{1};
    elseif iscell(Spks)
        for j = 1:length(Spks)
            if isempty(Spks{j})
            else
                if size(Spks{j}.times,1) == 1
                    Spks{j}.times = Spks{j}.times';
                end
                if ~isfield(Spks{j},'errs')
                    Spks{j}.errs = {};
                end
                AllData.AllClusters(j) = Spks{j};
                for k = 1:nspkt
                    if isfield(SpkTimes,'FullVData') && ~isempty(SpkTimes(k).FullVData)
                        AllData.FullVData = SpkTimes(k).FullVData;
                    end
                end
            end
        end
        
        if isfield(SpkTimes,'loadname')
            AllData.loadname = SpkTimes(1).loadname;
        end
        AllData.datenum = clusterdate;
        AllData.quickload = quickload;
    else
        AllData.Spikes = Spks;
        if ~isempty(AllSpks)
            AllData.AllSpks = AllSpks;
        end
    end
    if nspkt
        Expt.ClusterLoadTimes = [SpkTimes.loadtime];
    end
    Raw.stimlvl = stimlvl;
    Raw.stimch = stimch;
    if isstruct(UstimV)
        AllData.UstimV = UstimV;
    end
    if timeoffset
        Events.times = Events.times + timeoffset;
        f = {'Start' 'End' 'estimes' 'bstimes' 'TrueEnd' 'FalseStart'}
        for j = 1:length(f)
            if isfield(Expt.Trials,f{j})
                Expt.Trials.(f{j}) = Expt.Trials.(f{j})+timeoffset;
            end
        end
    end
       
    AllData.Events = Events;
    Expt.Probes = probes;
    if isfield(probes,'probe')
        Expt.Header.nprobes = max([probes.probe]);
    end
    if Expt.Header.Spike2Version < 1.23
        Expt = FixExpt(Expt,'ed');
    end
    iExpts = Expts;
    if isfield(Expt.Trials,'ve') && length(Expt.Trials.ve) < length(Expt.Trials.End)
        ve = mean(Expt.Trials.ve);
        Expt.Trials.ve(end:length(Expt.Trials.End)) = ve;
    end
    Expt.Trials = AddStimValsFromFile(strrep(name,'.mat','.SetStim'),Expt.Trials);
    if isfield(Expt.Trials,'badexpts')
        Expt.badexpt = Expt.Trials.badexpts;
        Expt.Trials = rmfield(Expt.Trials,'badexpts');
    end
    if isfield(Expts,'firsttrial')
            if timeoffset
        argon = {argon{:} 'timeoffset' timeoffset};
            end
            
            if ~isfield(Expts,'e3')
                [Expts.e3] = deal('e0');
            end
   [Expts, Expt] = SortExpts(Expts, Expt.Trials, Expt.Header, thecluster, Expt, state, argon{:});
    if ~exist('Exptlist','var')
        ExptList = MkExList(Expts);
        newlist = 0;
    else
        newlist = 1;
    end
   if Expt.newerrs || newlist  %% why not separate ifs? why rewrite idxfile, when newerrs?
       WriteErrors(idxfile, Expt);
       tExpts = Expts;
       Expts = iExpts;
       if timeoffset == 0 %don't re-write the idxfile
       save(idxfile,'Expt','Expts','ExptList');
       end
       Expts = tExpts;
   end
    end
    if isempty(Expts)
        nt = length(Expt.Trials.Start);
        [a,b] = Counts(Expt.Trials.Result);
        str = sprintf('No Expts in %s. %d Trials. Results',GetName(Expt.Header,'withsuffix'),nt); 
        for j = 1:length(a)
            str = sprintf('%s %d:%d ',str,b(j),a(j));
        end
        cprintf('red','%s\n',str);
    end
   if logfid >= 0
       fclose(logfid);
   end
   Expt.state = state;
   if saveexpts
       SaveExpts(name, Expts);
   end
    return;
end
state.resort = 1;  %force resort of expts if relisting
Expt.state = state;

%
% only get here if its the first time or called with 'relist'
if 0 & idxfile & exist(idxfile,'file')
    load(idxfile);
    return;
end

    if ~isfield(Expt,'errs')
        Expt.errs = errs;
        Expt.newerrs = length(errs);
    end
    if iscell(Spks)  && isfield(Spks{1},'values')
        AllData.AllSpikes = Spks;
        Spks = AllData.AllSpikes{1};
    elseif iscell(Spks)
        for j = 1:length(Spks)
            if isempty(Spks{j})
            else
                if size(Spks{j}.times,1) == 1
                    Spks{j}.times = Spks{j}.times';
                end
                AllData.AllClusters(j) = Spks{j};
            end
        end
        Spks = Spks{1};
    else
    AllData.Spikes = Spks;
        if ~isempty(AllSpks)
            AllData.AllSpks = AllSpks;
        end
    end



nt = 0;
nx = 0;
Expts = [];

%id = find(ismember(Spks.codes(:,1),thecluster));
Expt.setprobe = setprobe; 
frametimes = [];
bstimes = [];
estimes = [];

if exist(framechname,'var')
%On the PRoPix, we get up for 1 frame, down for the next
% so every event 
    if isfield(framech,'level')
%        id = find(framech.level == 1);
%        frametimes = framech.times(id) .* 10000;
         frametimes = framech.times .*10000;
    else
    frametimes = eval([framechname '.times * 10000']);
    end
    Header.frameperiod = median(diff(frametimes));
else
    Header.frameperiod = 167;
end

Header.AplayVer = state.version;
if ~isempty(ustimmarkch)
    ustimmarkch.times = ustimmarkch.times .* 10000;
else
    ustimes = [];
end

if isempty(stimlvl.level)
    cprintf('errors','Empty StimOn/Off\n');
    return;
end
fixlast = 0;
if fixup(1) && stimlvl.level(end) == 1 && stimlvl.level(1) == 0
    stimlvl.times = stimlvl.times(1:end-1);
    stimlvl.level = stimlvl.level(1:end-1);
    cprintf('errors','Deleting Last Stimlevel as it is ON\n');
    fixlast = 1;
end

if exist(stimch,'var') & isfield(stimlvl,'level') & length(stimlvl.times) > 1
%    id = find(Ch8.level == 1);
% need to add an extra event at the end of each list to avoid issues with
% empty finds on the last trial. But make it long after to avoid any
% confusion with the real one;
    maxt = max([max(Text.times) max(Events.times) max(stimlvl.times)]);

    if isfield(stimlvl,'inverted') && stimlvl.inverted
        bstimes = stimlvl.times(stimlvl.level == 0) * 10000;
        bstimes(end+1) = maxt+50000;
        estimes = stimlvl.times(stimlvl.level == 1) * 10000;
        estimes(end+1) = maxt+50010;
    else
        bstimes = stimlvl.times(stimlvl.level == 1) * 10000;
        bstimes(end+1) = maxt+50000;
        estimes = stimlvl.times(stimlvl.level == 0) * 10000;
        estimes(end+1) = maxt+50010;
    end
end


tic;

%use mat2cell NOT cellstr; Cellstr deblanks
aText.text =mat2cell(Text.text(:,1:end-1),ones(1,size(Text.text,1)),size(Text.text,2)-1);
id = find(strncmp('xxx=',aText.text,4));
ids = setdiff(1:size(Text.text,1),id);
aText.text = aText.text(ids);
xid = ones(size(id));
nl = size(Text.text,1)-length(id);
k=1;
m = 1;
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
        xid(m) = j;
        m = m+1;
% if we ever go back to this, probably need deblank(Text.....
    end
    if length(aText.text{k}) == 0 
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
    end
    k = k+1;
end

aText.times = Text.times(ids);
aText.codes = Text.codes(ids,:);
aText.text = deblank(aText.text);
fprintf('Reading Text took %.2f\n',toc);


[rfstr, rfdat] = MkUfl(name, aText);

Header.rfstr = rfstr;

bsid = find(strncmp('bs',aText.text,2));
bsstimes = aText.times(bsid);
esid = find(strncmp('es',aText.text,2));
esstimes = aText.times(esid);

tstore = zeros(size(bsid));
tstore(find(strncmp('bss',aText.text(bsid),3)))=1;
Trials.bsstimes = bsstimes;
Trials.esstimes = zeros(size(bsstimes));
for j = 1:length(esstimes)
    id = find(bsstimes < esstimes(j));
    if ~isempty(id)
    Trials.esstimes(id(end)) = esstimes(j);
    end
end
%Trials.esstimes(id) = NaN;

tic;
xid = find(strncmp('op=',aText.text,3)); %Text descriptors
opid = setdiff(find(strncmp('op',aText.text,2)),xid);

for j = 1:length(opid)
    str= sscanf(aText.text{opid(j)},'op%d');
    if ~isempty(str)
        storing(j) = str;
    else
        storing(j) = 0;
    end
end
storing = bitand(storing,STOREBIT);
%storeonoff is a list of times where storing is toggled
storeonoff = 1+find(abs(diff(storing)) > 0); 
if isempty(storeonoff)
    storeonoff = [1 1];
%what if storing was never on? e.g. for online data
    storetimes = [0 max(aText.times)];
else
    storetimes = [];
end

ve = find(strncmp('BGCS Version',aText.text,10));
if ve
    version = sscanf(aText.text{ve(1)},'BGCS Version %f');
else
    version = 1.1;
end

ids = find(strncmp('fz',aText.text,2));
for j = 1:length(ids)
fzs(j) = sscanf(aText.text{ids(j)},'fz%d');
end
if j > 0
    framerate = mean(fzs);
%    framerate = sscanf(Text.text(ids(1),:),'fz%n');
else
    if isfield(defaults,'fz')
        framerate = defaults.fz;
    else
        framerate = 96;
    end
end


%only read nfxx made by binoc internally, not nf=xxx by user
ids = find(strncmp('nf',aText.text,2) & ~strncmp('nf=',aText.text,3));
for j = 1:length(ids)
    fzs(j) = sscanf(aText.text{ids(j)},'nf%d');
end
if j > 0
    nomdur = prctile(fzs,90) .* 10000/framerate;
end

instim = 0;
inexpt = 0;
nextonoff = 1;
if onlinedata
    storestate = 1;
else
    storestate = storing(1);
    tonoff = Text.times(opid(storeonoff(nextonoff)));
end

    tic;
    Events.store = zeros(size(Events.times));
    if ~storing(storeonoff(1)) & storing(1) %% first event is an off
        onid = find(Events.times < Text.times(opid(storeonoff(1))));
        Events.store(onid) = 1;
    end
    
    for j = 1:length(storeonoff)
        if storing(storeonoff(j))
            if length(storeonoff) > j
                onid = find(Events.times >= Text.times(opid(storeonoff(j))) ...
                    & Events.times < Text.times(opid(storeonoff(j+1))));
            else
                onid = find(Events.times >= Text.times(opid(storeonoff(j))));
            end
            Events.store(onid) = 1;
        end
    end
    fprintf('Store Index takes %.2f\n',toc);
if onlinedata && sum(Events.store) == 0
    Events.store = ones(size(Events.times));
end

if onlinedata && sum(Events.codes(:,1)==STARTEXPT) == 0
    fprintf('Missing Expt Start First code at %.1f  ',Events.times(1)./10000);
    if logfid > 0
        fprintf(logfid,'Missing Expt Start First code at %.1f  ',Events.times(1)./10000);
    end
    id = find(Text.codes(:,1) == 1 & Text.times > 0.1);
    if ~isempty(id)
        Events.times = [Text.times(id(1)); Events.times];
        Events.store = [1; Events.store];
        Events.codes = [[STARTEXPT 0 0 0 ]; Events.codes];
        fprintf('True Start at %.1f\n',Text.times(id(1))./10000);
    end
end
nonstore = 0;
lastend = 0;
ntrials = sum(Events.codes(:,1) == FRAMESIGNAL);
trynew = 0;
if trynew %didn't help
Trials.Start(1:ntrials) = 0;
Trials.End(1:ntrials) = 0;
Trials.Trial(1:ntrials) = 0;
Trials.TrueEnd(1:ntrials) = 0;
Trials.Startev(1:ntrials) = 0;
Trials.stored(1:ntrials) = 0;
Trials.Result(1:ntrials) = 0;
Trials.serdelay(1:ntrials) = 0;
Trials.bstimes(1:ntrials) = 0;
Trials.delay(1:ntrials) = 0;
Trials.endelay(1:ntrials) = 0;
Trials.estimes(1:ntrials) = 0;
Trials.id(1:ntrials) = 0;
Trials.FalseStart(1:ntrials) = 0;
end

%estimes = estimes(estimes > 0);
%bstimes = bstimes(estimes > 0);
%remove stimON going to off at start of file 
if estimes(1) <= bstimes(1) && length(estimes) == length(bstimes)+1
    estimes = estimes(2:end);
elseif length(bstimes) == length(estimes)+1 %crash can leave trailing start e.g. ruf1914
    bstimes = bstimes(1:end-1);
end

txt = scanlines(strrep(name,'.mat','Add.txt'),'silent');
badexpts = [];
for j = 1:length(txt)
    if strfind(txt{j},'badexpt');
    badexpts(j) = sscanf(txt{id(j)},'%d');
    end
end

if ~isempty(stimchange)
    ds = diff(stimchange.times);
    [minds, id] = min(ds);
    if minds < 0.02;
        err = sprintf('Very Short Stimchange intervals at %.2f',stimchange.times(id));
        Expt = AddError( Expt, err);
    end
end

isi = bstimes(2:end)-estimes(1:end-1);

xid = find(Events.codes(:,1) == 255 & Events.codes(:,3) == 1); %delete SimOn
for j = 1:length(xid)
    ti = find(bstimes > Events.times(xid(j)),1);
    bstimes(ti) = [];
    ti = find(estimes > Events.times(xid(j)),1);
    estimes(ti) = [];
end

if length(bstimes) > 2
xid = find(Events.codes(:,1) ==  14 & Events.times > bstimes(1)+5000); %Trigger expt will make DIO pulse
for j = 1:length(xid)
    fprintf('Removing Pulse placed by TriggerExpt at %.3f\n',Events.times(xid(j)));
    [a,b] = min(abs(bstimes - Events.times(xid(j))));
    bstimes(b) = [];
    estimes(b) = [];
end
end

id = find(isi < 400); %shouldnt happen, but see lemM209.5
if length(id)
    nerr = nerr+1;
    err = sprintf('Bad ISIs (%d), starting at %.0f)\n',length(id),bstimes(id(1)));
    Expt = AddError(Expt, err);
    Result.badisi = id;
    fprintf('%d impossible isis\n',length(id));
% if this is a rapid pulse sequence, sometimes sent by binoc.
% then need to remove the first bstime preceding the first short isi
% N.B. this only works for first train
    if bstimes(id(1)+1)-bstimes(id(1)) < 200
        id = cat(1,id(1)-1,id);
    end
    useid = setdiff(1:length(bstimes), id+1);
    bstimes = bstimes(useid);
    estimes = estimes(useid);
    isi = bstimes(2:end)-estimes(1:end-1);
else
    Result.badisi = [];
end

xname = strrep(name,'.mat','Add.txt');
fid = fopen(xname,'r');
if fid > 0
    a = textscan(fid,'%f %s','delimiter','\n');
    fclose(fid);
    eid = find(strcmp('EndExpt',a{2}));
    for j = 1:length(eid)
        Events = Spike2.AddEvent(Events, a{1}(j), [ENDEXPT 0 0 4]);
    end
    sid = find(strcmp('AddStartExpt',a{2}));
    for j = 1:length(sid)
        Events = Spike2.AddEvent(Events, a{1}(j), [STARTEXPT 0 0 4]);
    end
    sid = find(strcmp('CombineExpts',a{2}));
    if ~isempty(sid)
        cprintf('red','%s Says to combine Expts\n',xname);
        state.combineexpts = 1;
    end
end

readmethod = 0;
settrials = 0;
if state.method == 1
    badbad = [];
    bnid = find(strncmp('bsn',aText.text,3));
    if onlinedata
        id = find(strncmp('bss',aText.text,3));
        if length(id) < length(bsstimes)/50 && length(id) < 5 %% online and storage off
            id = find(strncmp('bs',aText.text,2));
        end
    else
        id = find(strncmp('bss',aText.text,3));
    end
bid = id;
bssid = bid; 
bsstimes = aText.times(id);
    if onlinedata
        id = find(strncmp('es',aText.text,2));
    else
        id = find(strncmp('ess',aText.text,3));
    end
esstimes = aText.times(id); 
%Only find expts where storage was on at expt start. But any endexpt event
%will end, even if storage is off. E.G. a binoc crash. when binoc restarts
%can send EndExpt before setting storage on
exendid = find((Events.codes(:,1) == ENDEXPT | Events.codes(:,1) == CANCELEXPT));
et = Events.times(exendid);
allstartid = find(Events.codes(:,1) == STARTEXPT);
goodend = find(Events.codes(:,1) == ENDEXPT);
%When files are split, a slow close/open in Spike2 can lead to losing the
%START_EXPT CODE.  Detect by (1) more END thatn START and (2) missing start
%is early and (3) storage was on
if state.bysuffix && length(exendid) > length(allstartid)
    for j = 1:length(allstartid)
        id = find(exendid < allstartid(j));
%dont add back in a CANCELED Expt, or and EndExpt that is suspiciously soon        
        if length(id) >= j  && exendid(id(end)) > 100 && Events.codes(exendid(id(end)),1) == ENDEXPT
            if id(end) > 1
                xid = exendid(id(end-1)):exendid(id(end)); %Event codes between two ENDEXPT
                if sum(Events.store(xid)) > length(xid) * 0.8 %storage was on
                    Events = Spike2.AddEvent(Events, Events.times(xid(2))-0.1, STARTEXPT)
                    allstartid = cat(1,allstartid(1:j-1), xid(2), allstartid(j:end));
                    Expt = AddError( Expt, sprintf('Missing StartExpt at %.3f. Replacing',Events.times(xid(2))./10000));
                    state.combineexpts = 2;
                end
            elseif Events.codes(exendid(1),1) == ENDEXPT %missing start at the beggining. but have end then  start then end later.
                xid = find(Events.codes(:,1) == STARTSTIM);
                Events = Spike2.AddEvent(Events, Events.times(xid(1))-0.5, STARTEXPT)
                allstartid = cat(1, xid(1), allstartid);
                Expt = AddError( Expt, sprintf('Missing StartExpt near %.3f. Replacing',Events.times(xid(2))./10000));
                state.combineexpts = 2;
            end
        end
    end
    
end
for j = 1:length(allstartid)
    if j < length(allstartid)
    id = find(et > Events.times(allstartid(j) ) & et < Events.times(allstartid(j+1)));
    else
    id = find(et > Events.times(allstartid(j)));
    end
    if ~isempty(id)
        eid = exendid(id(1));
        nsaved = sum(bsstimes > Events.times(allstartid(j)) & bsstimes < Events.times(eid));
        if nsaved > 2 && Events.store(allstartid(j)) == 0
           cprintf('blue','Storage Turned on Mid expt %.2f-%.2f\n',Events.times(allstartid),Events.times(eid));
            Events.store(allstartid(j)) = 1;
        end
    end
end
unsavedexpts = sum(Events.codes(:,1) == STARTEXPT & Events.store == 0);
exstartid = find(Events.codes(:,1) == STARTEXPT & Events.store > 0);
exendid = find((Events.codes(:,1) == ENDEXPT | Events.codes(:,1) == CANCELEXPT));
%exendid = find((Events.codes(:,1) == ENDEXPT | Events.codes(:,1) == CANCELEXPT) & Events.store > 0);
texstartid = find(Text.codes(:,1) == STARTEXPT);

%bstimes has an artificial extra on the end
if length(bsstimes) < length(bstimes)-1
    t = bsstimes(1)-bstimes;
%only 1 event from -1000ms to + 100ms) from text mark = safe to assume any before this are related to file open
    id = find(t > -1000 & t  <10000); 
    if isempty(id)
        Expt.t = bsstimes(1);
        Expt = AddError(Expt, 'First Stimon %.2f before first bss\n',t(1));
        if t(1) > 100000
            id = find(t > 10000);
            bstimes = bstimes(id(end)+1:end);
            estimes = estimes(id(end)+1:end);
            Expt = AddError(Expt, '%s Stimons first bss',length(id));
        end            
    else
        Expt.t = bsstimes(1);
        if (id(end) > 3)
            Expt = AddError(Expt, '%d StimOns before first saved trial (%.2f)\n',id(end),bsstimes(1));
        else
            cprintf('blue', '%d StimOns before first saved trial (%.2f)\n',id(end));
        end
        bstimes = bstimes(id:end);
        estimes = estimes(id:end);
    end
end
evid = [];
bsid = [];
allfsid = [];
ExptStart = Events.times(exstartid);
exendid = exendid(exendid > exstartid(1));
ExptEnd = Events.times(exendid);
if isempty(ExptEnd) %% sometimes .mat file is missing EndExpt Marker. But best to fix .smr
    Expt = AddError(Expt, 'No EndExpt Marker Expt');
    ExptEnd = Events.times(end);
    exendid = length(Events.times)+1;
    Events.codes(exendid,1) = ENDEXPT;
    Events.store(exendid) = Events.store(end);
    Events.times(exendid) = Events.times(end)+0.1;
end
ExptCode = Events.codes(exendid,1);
ExptStart(length(exstartid)+1) = ExptEnd(end) +1;

for j = 1:length(exstartid)
    ts = ExptStart(j);
    id = find(ExptEnd > ExptStart(j)  & ExptEnd < ExptStart(j+1));
    if isempty(id) %force end Expt if notthing - probabaly a crash
        if j < length(exstartid)
            te = ExptStart(j+1)-10000;
            Expts(j).result = ENDEXPT;
        else
        te = Events.times(end);
        Expts(j).result = ENDEXPT;
        end
    else
        te = ExptEnd(id(1));
        Expts(j).result = ExptCode(id(1));
    end
    id = find(Events.times > ts & Events.times < te);
    evid = cat(1,evid,id);
%take care of any hanging estime after expt end. Can happen with
%fixed/crashed files
    eid = find(estimes > ts & estimes < te);
    if isempty(eid)  || ismember(j,badexpts) %no trials in this expt
        Expts(j).result = -1;
        Expts(j).storecount = 0;
    else
        %if final trial in expt is a badfix, and its short, then final FRAMESIGNAL
        %code in serial line can be AFTER the StimOn Channel goes to off.  So
        %search for all FRAMESIGNALS before end Expt (no more esttimes anyway
        %becuase of previous test.
        fsid = find(Events.codes(id,1) ==  5 & Events.times(id) < te);
        fstimes = Events.times(id(fsid));
        bid = find(bstimes > ts & bstimes < te);
        %
        %fstimes is times of FRAME_SIGNAL codes. fsid thier index in Events
        %bstimes(bid) is StimChan On events in expt
        %estimes(eid) is StimChan Off events in expt
        
        bids{j} = bid;
        fsids{j} = id(fsid);
        bsid = cat(1, bsid, bid);
        allfsid = cat(1, allfsid, id(fsid));
        Result.startcounts = [length(fsid) length(bid)];
        Result.storecount = sum(Events.store(fsid));
        Result.name = name;
        if ~isempty(fsid) && ~isempty(bid)
            td = (bstimes(bid(end))-Events.times(id(fsid(end))))./10000;
        end
        if sum(Result.startcounts == 0) %either fsid or bid is empty
            fprintf('No Trials for block %d %.1f - %.1f\n',j,ts,te);
            Result.bsdelay = [];
            fsids{j} = [];
            bids{j} = [];
            bsid = bsid(1:end-length(bid));
        elseif length(fsid)+1 == length(bid) &&  td > 10
            %if there is StimON at the end of hte expt that trails by a long way, probaly an error
            %E.G. jbeG016 at 7496.0
            cprintf('blue','Ex %d Unmatched final Stimlevel at %.1f, %.1f after last bss\n',j,bstimes(bid(end))./10000,td);
            bid = bid(1:end-1);
            bsidx(id(fsid)) = bid;
            Result.bsdelay = Events.times(id(fsid)) - bstimes(bid);
            bids{j} = bid;
            bsid = bsid(1:end-1);
        elseif length(fsid) == length(bid)
            fprintf('Ex %d Stimlevel and FrameSignal Lengths match %d(%.1f-%.1f)\n',j,length(fsid),bstimes(bid(1)),bstimes(bid(end)));
            bsidx(id(fsid)) = bid;
            Result.bsdelay = Events.times(id(fsid)) - bstimes(bid);
        elseif length(fsid) == length(bid)+1 && fixlast && fstimes(end) > estimes(bid(end)) %extra frame signal belongs to chopped final StimOn
            fprintf('Ex %d Removing Trailing Frame signal %d(%.1f)\n',j,length(fsid),Events.times(id(fsid(end))));
            fsid = fsid(1:end-1);
            allfsid = allfsid(1:end-1);
            fsids{j} = fsid;
            bids{j} = bid;
            Result.bsdelay = Events.times(id(fsid)) - bstimes(bid);
        else
            fprintf('Length Mismatch for Stimlevel (%d) and FrameSignal (%d) Last gap %.1f\n',length(bid),length(fsid),td);
            [a,b] = fileparts(name);
            GetFigure(b);
            s = FindMissingTimes(Events, aText, bstimes, estimes,bid,id(fsid),ts,te);
            if ~isempty(s)
                Expt = AddError(Expt,s);
            end
            Result.bsdelay = [];
        end
        Expts(j).storecount = sum(Events.store(fsids{j}));
    end
    
    Expts(j).start = ts;
    Expts(j).end = te;
    if length(bid) && Expts(j).result >= 0
        Expts(j).starti = bid(1);
        Expts(j).endi = bid(end);
    else
        Expts(j).starti = NaN;
        Expts(j).endi = NaN;
    end
end
%Expts = Expts([Expts.result] ==ENDEXPT);
evid = unique(evid);
bsid = unique(bsid);
fsid = find(Events.codes(evid,1) ==  5); %frame signal
fstimes = Events.times(evid(fsid));
fsid = allfsid;
fstimes = Events.times(fsid);
for j = 1:length(Expts)
    Expts(j).firsttrial = find(bsid == Expts(j).starti);
    Expts(j).lasttrial = find(bsid == Expts(j).endi);
    if isempty(Expts(j).lasttrial) && ~isempty(Expts(j).firsttrial)
        cprintf('red','Empty Last Trial Expt %d\n',j);
    end
    good(j) =  ~isempty((Expts(j).firsttrial));
end
rwid = find(Events.codes(:,1) == 'q');
rwtimes = Events.times(rwid);
Expts = Expts(good);

if isempty(Expts)
    fprintf('No Expts in %s\n',name);
    return;
end
bfid = find(aText.codes(:,1) == BADFIX & aText.codes(:,4) ==2);
offid = find(aText.codes(:,1) == AUTOSTOREOFF & aText.codes(:,4) ==2);
bfevid = find(Events.codes(:,1) == BADFIX);
endtrid = find(Events.codes(:,1) == ENDTRIAL);
sacid = find(strncmp('Sa:',aText.text,3)); %lines sent by binoc if ending trial
endstimid = find(strncmp('EndStim',aText.text,7));
t = mygetCurrentTask();
if t.ID == 0
    guiok = 1;
else
    guiok = 0;
end





Trials.Start = zeros(size(fsid))';
Trials.End = zeros(size(fsid))';
Trials.FalseStart = zeros(size(fsid))';
Trials.TrueEnd = zeros(size(fsid))' .* NaN;
Trials.delay = Trials.TrueEnd;
ts = now;
if state.profiling == 2
    profile on;
end
if length(fsid) == length(bsid)
   Result.bsdelay = fstimes - bstimes(bsid);
   readmethod = 1;
   if guiok
       wh = waitbar(0,'Building Trial Time List');
   end
%Trials.EndTxt idetnifies the text event that matched the end of trial n   
  for j = 1:length(fsid)
      Trials.Start(j) = bstimes(bsid(j));
      Trials.End(j) = estimes(bsid(j));
      Trials.bstimes(j) = Trials.Start(j);
      Trials.estimes(j) = Trials.End(j);
      Trials.Trial(j) = j+starttrial;
      Trials.Result(j) = 1;
      Trials.bsdelay(j) = fstimes(j)-bstimes(bsid(j));
      id = find(esstimes > fstimes(j),1);
      if length(id) == 0 %% all remaining trials are bad fix, will get here
          Trials.Result(j) = -1;
         Trials.EndTxt(j) = 0;
     elseif j <= length(fsid)
         Trials.EndTxt(j) = esstimes(id(1));
         if j < length(fsid) && esstimes(id(1)) > fstimes(j+1)
             a = find(aText.times(bfid) > fstimes(j) & aText.times(bfid) < fstimes(j+1));
             if ~isempty(a)  %Bad fixation trial
                 Trials.Result(j) = 0;
                 Trials.EndTxt(j) = aText.times(bfid(a(1)));
             else
                 a = find(aText.times(endstimid) > fstimes(j) & aText.times(endstimid) < fstimes(j+1));
                 if ~isempty(a)  %Also Bad fixation trial
                     Trials.Result(j) = 0;
                     Trials.EndTxt(j) = aText.times(endstimid(a(1)));
                 else
                     a = find(Events.times(bfevid) < fstimes(j+1));
                     Trials.Result(j) = -1;
                     if ~isempty(a) && Events.times(bfevid(a(end))) > fstimes(j) %Bad fixation trial
                         bftime = Events.times(bfevid(a(end)));
                         a = find(aText.times > bftime,1);
                         Trials.EndTxt(j) = aText.times(a(1));
                         Trials.Result(j) = -2;
                         a = find(aText.times(sacid) > fstimes(j) & aText.times(sacid) < fstimes(j+2));
                         if isempty(a)
                             fprintf('Badfix at %.0f, no text\n',bftime);
                         else
                             x = aText.text{sacid(a(1))};
                         end
                     else
                         a = find(Events.times(endtrid) < fstimes(j+1) & Events.times(endtrid) > fstimes(j));
                         if ~isempty(a)
                             bftime = Events.times(endtrid(a(1)));
                             a = find(aText.times > bftime);
                             Trials.EndTxt(j) = aText.times(a(1));
                             Trials.Result(j) = -3;
                             if state.alltrials
                                 fprintf('Trial end no text %.0f\n',bftime);
                             end
                         else %a mess no end trial between starts j and j+1
                             a = find(aText.times(offid) < fstimes(j+1) & aText.times(offid)  > fstimes(j));
                             if ~isempty(a)
                                 bftime = aText.times(offid(a(1)));
                                 Trials.EndTxt(j) = bftime;
                                 Trials.Result(j) = -4;
                                 fprintf('Sorage off after Badfix without report %.0f\n',bftime);
                             end
                         end
                     end
                 end
             end
         end
         else
         Trials.EndTxt(j) = 0;
      end
      if j > 1 && Trials.EndTxt(j) < Trials.EndTxt(j-1) && Trials.EndTxt(j)> 0
          x = Trials.EndTxt(j) - Trials.EndTxt(j-1);
%when this happens, it can serioulsy mess up alignment of trial data later
          cprintf('red','End Trial text timing error at %.0f. SERIOUS!\n',Trials.EndTxt(j));
      end

      if j < length(fstimes)
          nextfs = fstimes(j+1);
      else
          nextfs = estimes(end)+2;
      end
      id = find(rwtimes > Trials.Start(j) & rwtimes < nextfs);
      if ~isempty(id)
          Trials.rwtimes{j} = rwtimes(id);
      else
          if Trials.estimes(j)-Trials.bstimes(j) > 20000 & Trials.Result(j) > 0
              Trials.rwtimes{j} = 0;
          else
              Trials.rwtimes{j} = 0;
          end
      end
      if guiok && mod(j,1000) == 0
          waitbar(j/length(fstimes));
      end
  end
  if guiok
      delete(waitbar(1));
      st = now;
      wh = waitbar(0,'Matching Frame Times');
      drawnow;
      tt = mytoc(ts);
  end
  if ~isempty(frametimes) % have VTR channel
      
      ff = 1;
      nf = length(frametimes);
      maxl(2:length(Trials.Start)) = 100+ceil(diff(Trials.Start)./100);
      maxl(1) = ceil(Trials.Start(1)./100);
%usnig bsxfun is slower.        
%      delays = bsxfun(@minus, frametimes, Trials.Start);
      ix = 1;
      Start = Trials.Start;
      Ends = Trials.End;
      Trials.endelay = zeros(size(Trials.End));
      for j = 1:length(Start)
%              delays = frametimes-Start(j);
              id = find(frametimes- Start(j) > 0,1);
%          id = find(delays(ix:end,j) > 0 & delays(ix:end,j) < 20000);
          if ~isempty(id) && frametimes(id(1)) - Start(j) < 20000
              Trials.delay(j) = frametimes(id(1)) - Trials.Start(j);
              Trials.Start(j) = frametimes(id(1));
              Trials.FalseStart(j) = 0;
              ff = id(1);
              last = min([ff+1000 nf]);
              ix = id(1);
          else
              Trials.FalseStart(j) = 1;
              Trials.delay(j) = NaN;
              last = nf;
          end
          id = find(frametimes > Trials.End(j),1);
          if ~isempty(id) & frametimes(id(1))-Trials.End(j) < 500
              Trials.endelay(j) = Trials.End(j) - frametimes(id(1));
              Trials.End(j) = frametimes(id(1));
              Trials.TrueEnd(j) = frametimes(id(1));
          else
              fprintf('No Frametimes after End of Trial at %.2f\n',Ends(j));
          end
          if guiok && mod(j,1000) == 0
              waitbar(j/length(fstimes));
          end
      end
  end
  if guiok
      delete(waitbar(1));
      drawnow;
      fprintf('Frames took %.2f\n',mytoc(ts)-tt);
  end
  %badfix is always (?? what about Sa?) detected by spike2, so the stimoff
  %marker will always be after this.
  bid = find(Events.codes(:,1) == 11); %bad fix
  for j = 1:length(bid)
      id = find(Trials.bstimes < Events.times(bid(j)));
      if length(id)
      Trials.Result(id(end)) = 0;
      bfdelay(id(end)) = Trials.estimes(id(end))-Events.times(bid(j));
      else
          badbad = [badbad Events.times(bid(j))];
      end
  end
  bid = find(aText.codes(:,1) == 11 & aText.codes(:,4) == 2); %bad fix from spike2
  gbid = find(strncmp('BAD Saccade',aText.text(bid),10));  %Did complete trial, but invalid psych sacc
  bid = setdiff(bid, bid(gbid));
  for j = 1:length(bid)
      id = find(Trials.bstimes < aText.times(bid(j)));
      if length(id)
          bfdelay(id(end)) = Trials.estimes(id(end))-aText.times(bid(j));
          if bfdelay(id(end)) < 0 && Trials.Result(id(end)) > 0 %did complete before badfix
              fprintf('Trial %d Badfix after Stimulus completed at %.3f\n',id(end),Trials.estimes(id(end))./10000);
              Trials.Result(id(end)) = 2;
          else
              Trials.Result(id(end)) = 0;
          end
      else
          badbad = [badbad aText.times(bid(j))];
      end
  end
      


  settrials = 1;
  nt = length(fsid);
else
    badcount = [];
    for j =1:length(fsids)
        if length(fsids{j}) ~= length(bids{j})
            badcount(j) = j;
        end
    end
    Expt = AddError(Expt,'Cant use New Read method - mismatched counts at%s\n',sprintf(' %d',badcount));
    readmethod = -1;
end
end

tt = mytoc(ts);




if settrials == 0 %were not set by new method
    fprintf('Using old Read Method\n');
    readmethod = 0;
    id = find(Text.codes(:,1) == BADFIX);
    bftimes = Text.times(id);
    for j = 1:size(Events.codes,1)

% if storage turned of mid-stim, don't want to miss ENSTIM marker
    if Events.codes(j,1) ~= ENDSTIM
        storestate = Events.store(j);
    end        
    if Events.codes(j,1) == FRAMESIGNAL
        nowt = Events.times(j);
        if storestate
            if nt & Trials.Result(nt) < 0
                nt = nt;
            end
            nt = nt+1;
        Trials.Start(nt) = Events.times(j);
        Trials.End(nt) = Trials.Start(nt)+nomdur; %% just in case an online file is missing end
        if findtrial  & Trials.Start(nt) > findtrial
            findtrial = 0;
        end
        Trials.Startev(nt) = Events.times(j);
        Trials.stored(nt) = storestate;
        Trials.Result(nt) = 1;
        if nt > 1 & length(Trials.End) < nt-1
            nerr = nerr+1;
            err = sprintf('Missing end Trial %d (EX %.0f, start %.0f)\n',nt-1,inexpt,Trials.Start(nt-1));
            Expt = AddError(Expt, err);
        end
        Trials.Trial(nt) = nt+starttrial;
        instim = 1;
%the event time can be just before the stim signal since it is not delayed
%to the sync (? correct explanation. Definitely what happens
% is the times of teh StimulusON digital event markers
%frametimes are the Vsync digital markers (every frame)
% the digital step can be >20ms after the serial signal is received, e.g.
% in ruf1989 at 209.43 sec
% in lem017 at 120.6 it is nearl 80ms late. at 289.6 its 300ms late. Looks
% like we could check for immediately preceding STARTSTIM (6) being before 
% StimON to check for this

        id = find(bstimes < Events.times(j)+300);
        
        if Events.times(j) > 238170000
            bstimes(id(end));
        end
        if id
            Trials.serdelay(nt) = Events.times(j) - bstimes(id(end));
            Trials.bstimes(nt) = bstimes(id(end));
        end
        if ~isempty(id) & bstimes(id(end)) > lastend 
            if id(end) < length(bstimes) & bstimes(id(end)+1) - Trials.Start(nt) < 10 %< 1ms to next = probably early
                Trials.Start(nt) = bstimes(id(end)+1);
            elseif bstimes(id(end)) > lastend
                Trials.Start(nt) = bstimes(id(end));
            end
            if Trials.serdelay(nt) > 10000
                nt = nt;
            end
            Trials.stored(nt) = storestate;
%here, Trials.Start is the time of the Digital event marker
%If vertical retrace was recorded, set start time to next one of these
%To a first approximation, it is the find(framemtimes ......  calls (here
%and for ENDSTIM) that take all the time in this loop. 
            if ~isempty(frametimes) % have VTR channel
                id = find(frametimes > Trials.Start(nt) & frametimes < Trials.Start(nt)+500);
                if ~isempty(id)
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                    Trials.FalseStart(nt) = 0;
                else
                    Trials.FalseStart(nt) = 1;
                    Trials.delay(nt) = NaN;
                end
            else
%                Trials.Start(nt) = Events.times(j);
                Trials.FalseStart(nt) = 2;
                Trials.delay(nt) = NaN;
            end
        else
            Trials.Start(nt) = Events.times(j);
            Trials.delay(nt) = NaN;
            if isempty(id)
                Trials.FalseStart(nt) = 1;
            elseif id(end) < length(bstimes) & bstimes(id(end)+1) - Events.times(j) < 800 ...
%                    & Events.codes(j-1) == STARTSTIM ... %< 1ms to next = probably early
                err = sprintf('StimON at %.2f is %.1f ms late but STARTSTIM at %.2f',...
                bstimes(id(end)+1),(bstimes(id(end)+1)-Events.times(j))./10,Events.times(j-1));
                Expt = AddError(Expt, err);
                Trials.FalseStart(nt) = 0;
                
            else
%Serial input can be very late if Spike2 got busy. Use the DIO stimon -
%this is the true start
                Trials.FalseStart(nt) = Events.times(j) - bstimes(id(end));
                Trials.Start(nt) = bstimes(id(end));
%This error also generated if final trial doesn not ed (storage stopped
%first
                 if length(bstimes) - id(end) < 2
                     err = sprintf('Missing Final StimON at %.2f(ev) %.2f(bss), but STARTSTIM at %.2f',Events.times(j),bstimes(id(end)),Events.times(j-1));
                     Expt = AddError(Expt, '-silent', err);
                 else
                     err = sprintf('Missing StimON at %.2f(ev) %.2f(bss), but STARTSTIM at %.2f',Events.times(j),bstimes(id(end)),Events.times(j-1));
                     Expt = AddError(Expt, err);
                 end
                id = find(frametimes > Trials.Start(nt));
                if ~isempty(id) 
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                end
            end
            Trials.stored(nt) = storestate;
        end 
        else
            nonstore = nonstore+1;
        end %if storestate
  
    elseif Events.codes(j,1) == ENDSTIM & storestate & nt
    if abs(Events.times(j) - 62740381) < 200
            Trials.End(nt) = Events.times(j);
            Trials.endelay(nt) = NaN;
    end
        if estimes
            id = find(estimes < Events.times(j)+500);
            if(id)
                Trials.TrueEnd(nt) = estimes(id(end));
                Trials.End(nt) = estimes(id(end));
                Trials.endelay(nt) = NaN;
                Trials.estimes(nt) = estimes(id(end));
  %if this is out by 400ms, probably failed to find correct end mark
                if Trials.TrueEnd(nt) < Events.times(j) - 4000  
                    fprintf('End event %.3f but marker %.3f\n',...
                        Events.times(j)./10000,Trials.TrueEnd(nt)./10000);
                    if Trials.End(nt) < Trials.Start(nt) & length(estimes) > id(end) & ...
                        estimes(id(end)+1) - Events.times(j) < 10000
                        Trials.TrueEnd(nt) = estimes(id(end));
                        Trials.End(nt) = estimes(id(end));
                        Trials.endelay(nt) = NaN;
                        Trials.estimes(nt) = estimes(id(end));
                    end
                end

                if ~isempty(frametimes) % have VTR channel
                    id = find(frametimes > Trials.End(nt));
                    if ~isempty(id) & frametimes(id(1))-Trials.TrueEnd(nt) < 500
                        Trials.endelay(nt) = Trials.End(nt) - frametimes(id(1));
                        Trials.End(nt) = frametimes(id(1));
                        Trials.TrueEnd(nt) = frametimes(id(1));
                    else
                    end
                end
            else
                Trials.TrueEnd(nt) = 0;
            end
            
        end
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 1;
        if (Trials.End(nt) - Trials.Start(nt)) < 1000
            instim = 0;
        end
        instim = 0;
        if Trials.TrueEnd(nt)
            lastend = Trials.TrueEnd(nt);
        else
            lastend = Trials.End(nt);
        end
    elseif Events.codes(j,1) == ENDTRIAL & storestate
        if instim
%  can't figure this out here because the BADFIX is only recorded in text,
%  not SampleKey (becuase this is send from Spike2, not received by, and
%  setting codes for sample keys is such a pain. But maybe should make all
%  of these events with code2 set to indicate it is from Spike2?
% Seems like this happens when fixation is broken just BEFORE stimulus on,
% but Spike2 has not registered this yet e.g. ruf2000 at 8867.9
%So, check badfix times from text (bftimes). If there is a badfix before
%this event, but after Stim start, nothing to worry about.
           id = find(bftimes < Events.times(j) & bftimes > Trials.Start(nt));
           Trials.Result(nt) = -1;  % this will be set to 0 if a BadFix is found.
           if isempty(id) 
               bsid = find(bstimes < Events.times(j));
               esid = find(estimes < Events.times(j));
               if length(bsid) && length(esid)
                   dur = estimes(esid(end))-bstimes(bsid(end));
               else
                   dur = (Events.times(j) - Trials.Start(nt));
               end
               if dur > 1000
                   err = sprintf('End Trial without End stim: %d (%.2f) dur %.1fms',...
                       nt-1,Events.times(j)/10000,(Events.times(j) - Trials.Start(nt))/10,dur/10);
                   Expt = AddError(Expt, err);
                   Trials.Result(nt) = -2;  % this will be set to 0 if a BadFix is found.
               end
           end
            Trials.End(nt) = Events.times(j);
            Trials.TrueEnd(nt) = NaN;
        end
    elseif Events.codes(j,1) == BADFIX & storestate %% Doesn't happen. Badfix is in Text, because it is sent, not received
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 0;
        instim = 0;
    elseif Events.codes(j,1) == STARTEXPT & storestate
        if inexpt %close an existing expt (e.g. if crashed out)
            Expts(nx).end = Events.times(j);
            Expts(nx).lasttrial = nt;
        end
        nx = nx+1;
        Expts(nx).start = Events.times(j);
        Expts(nx).firsttrial = nt+1;
        inexpt = 1;
    elseif Events.codes(j,1) == STARTEXPT % non stored
        inexpt = 0;
    elseif Events.codes(j,1) == ENDEXPT & nx & inexpt %storage was on at start expt
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        inexpt = 0;
        Expts(nx).result = ENDEXPT;
    elseif Events.codes(j,1) == CANCELEXPT & nx && inexpt %Cancel cancels an expt, even if storage off
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        Expts(nx).result = CANCELEXPT;
        inexpt = 0;
    elseif Events.codes(j,1) == ENDEXPT
        nx = nx;        
    elseif Events.codes(j,1) == CANCELEXPT
        nx = nx;        
    end
end

end
Trials.id = zeros(size(Trials.Trial))'; %needs to be a row 
if nt == 0
    Expts = [];
    return;
end
if inexpt
    Expts(nx).lasttrial =nt;
end
if length(Expts) > nx  && settrials == 0 %some Expts set in first pass never found
    fprintf('Only found %d of %d good Expts\n',nx,length(Expts));
end
ntrials = nt;
fprintf('Setting Trials Took %.2f\n',mytoc(ts));
tic;
if state.profiling == 2
    return;
end

trial = 1;
k = 1;
ix = 1;
nx = 1;


nmiss = sum(Trials.FalseStart ==1);
if nmiss > sum(tstore ==0)
    Expt = AddError(Expt,'Missing Frame time for %d trials',nmiss);
end    
%kludge. Should fill esstimes with when the badfix happened.
if length(Trials.End) > length(Trials.esstimes)
    nmiss = length(Trials.End) - length(Trials.esstimes)
    Expt = AddError(Expt,'Removing %d trials with no stimlevel in trial',nmiss);
    if readmethod == 1
        Trials.End = Trials.End(1:length(Trials.esstimes));
    else        
        Trials.End = Trials.End(1:length(Trials.esstimes));
        Trials.Start = Trials.Start(1:length(Trials.esstimes));
    end
end
id = find(Trials.esstimes(1:length(Trials.End)) == 0);
Trials.esstimes(id) = Trials.End(id);

% trynew should work, but need to check with some xxx= strings.
%problem with new method was cell2str deblanked. mat2cell works better
%

trynew = 1;
if trynew
end

for j = 1:length(Expts)
    Expts(j).midtrial = (Expts(j).lasttrial+Expts(j).firsttrial)/2;
end

ix = 1;
if trynew == 0
tic;
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
    else
% if we ever go back to this, probably need deblank(Text.....
    aText.text{k} = Text.text(j,1:end-1);
    aText.times(k) = Text.times(j);
    aText.codes(k,:) = Text.codes(j,:);
    end
    if j == 206991
        k
    end
    if length(aText.text{k}) %%? safe may remove blank strings that have codes.
        k = k+1;
    else
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
        k = k+1;
    end
end
aText.text = deblank(aText.text);
fprintf('Reading Text %.2f\n',toc);
end
%
%AddTxtFile allows problems with data files to be fixed by adding
%additional lines written by hand. Fornmat is
%time  text
%where time is an int in timestamp units (0.1ms)

[aText, Text] = AddText(regexprep(name,'\.[0-9]*.mat','Add.txt'), aText, Text);
[aText, Text,newlines] = AddText(strrep(name,'.mat','Add.txt'), aText, Text);

%set codes(:,3) to indicate storate status if not set in spike2
if isempty(storetimes) && sum(aText.codes(:,3) ==0)
    if storing(storeonoff(1)) 
        ss = 1;
    else
        ss = 2;
    end
    if length(storeonoff) == 2
        storeid(1,:) = opid(storeonoff);
        storetimes(1,1) = aText.times(opid(storeonoff(ss)));
    elseif length(storeonoff) == 1        
        if storing(storeonoff) ==0 %only change is on ->off
            storeid(1,2) = opid(storeonoff);
            storeid(1,1) = 1;
        else
            storeid(1,1) = opid(storeonoff);
        end
        if ss == 1
            storetimes(1,1) = aText.times(opid(storeonoff(1)));
        end
    end
    for j = ss:2:length(storeonoff)-1
        ne = ceil(j/2);
        storetimes(ne,1) = aText.times(opid(storeonoff(j)));
        storeid(ne,1) = opid(storeonoff(j));
        storetimes(ne,2) = aText.times(opid(storeonoff(j+1)));
        storeid(ne,2) = opid(storeonoff(j+1));
    end
    if storing(end) > 0
        storeid(end,2) = length(aText.times);
        storetimes(end,2) = aText.times(end);
    end
    for j = 1:size(storeid,1)
        aText.codes(storeid(j,1):storeid(j,2),3) = 1;
    end
    Expt.storetimes = storetimes;
end



Trials.Stimseq = {};
Trials.op = 0;
intrial = 0;
Peninfo.trode = '';
nrw = 0;
tic;
tstart = now;
lasttook = 0;
Stimulus.CorLoop = 0;
Stimulus.SpikeGain = 50; %default
Stimulus.id = 0;
StimTypes.num.id = 1;
StimTypes.num.SpikeGain = 1;
StimTypes.num.CorLoop = 1;
Stimulus.vs=0;
Stimulus.sq=0;

gotend = 0;
txtid = [];
lastfix.fx = 0;
lastfix.fy = 0;
fprintf('Text->stims .....',toc);
if guiok
    waitbar(0,sprintf('parsing %d txt lines',length(aText.text)));
end
if state.profiling == 1
    profile on;
end


id = find(strncmp('exptlabel',aText.text,8));
for j = 1:length(id)
    aText.text{id(j)} = strrep(aText.text{id(j)},'explabel','exptlabel');
end
id = find(strncmp('expvars',aText.text,7));
for j = 1:length(id)
    aText.text{id(j)} = strrep(aText.text{id(j)},'expvars','exptvars');
end

id = [];
codes = zeros(1,length(aText.text));
findstrs = {'puA' 'puF' 'USd' 'USp' 'USf' 'nph' 'ijump' 'mixac' 'baddir' ...
    'e1max' 'backMov' 'FakeSig' 'pBlack' 'aOp' 'aPp' 'seof' 'serange' ...
    'nimplaces' 'usenewdirs' 'choicedur' 'cha' 'imi' 'choicedur' 'ePr' 'prev' ...
    'coarsemm' 'psyv' 'imi' 'nbars' 'imY' 'imX' 'bjv' 'imx'  'imy' 'bpos' 'Trw'};
%charstrs are character variables that might be needed trial by tial
charstrs = {'Covariate'  'hxtype' 'cx' 'adapter' 'exp' 'et' 'e2' 'e3' 'exptlabel' 'exptvars' 'imprefix' 'stimtag' 've' 'Bs' 'cmdfile' 'VisualArea' 'Vn'};
%extrastrs at variabes we won't typically want trial by trial
%So its faster to find them now, and test each line for membership
extrastrs = {'StartDepth' 'Electrode' 'Write' ' Electrode' 'Experiment' 'Off at' 'CLOOP' 'cm=' ...
    'mt' 'fx' 'fy' 'fp' 'rw' 'st' 'fl+' 'Nf' 'dx:' 'EndExpt' 'sonull' 'NewConnect' 'BGCS Version' 'rptframes' ...
    'expname' 'immode' 'rcfile' 'BadFix' 'useexvals' 'centerxy' ...
    'Tube Prot',...
    'seof' '/local' 'imve' 'Sa:' 'op' 'annTyp' 'uf' 'lo' 'vve' 'testflag' 'Hemisphere' 'Reopened' 'monitor' 'RightHemi' 'ui' 'xvals' 'ch10' 'ch11';};
%really need to order these ny length, logest last, so that long matches
%take preceendnce
vartypes(1:length(findstrs)) = 'N';
vartypes((1+length(findstrs)):(length(findstrs)+length(charstrs))) = 'C';
vartypes((1+length(vartypes)):(length(vartypes)+length(extrastrs))) = 'X';
findstrs = [findstrs charstrs extrastrs];
for j = 1:length(findstrs)
    lens(j) = length(findstrs{j});
end
[a,b] = sort(lens);
findstrs = findstrs(b);
vartypes = vartypes(b);
for j = 1:length(findstrs)
    f = findstrs{j};
    slens(j) = length(f);
    id = find(strncmp(f,aText.text,length(f)));
    codes(id) = j;
end

id = find(strcmp('uf',findstrs));
id = find(codes ==id);
if ~isempty(id)
    ufline = aText.text{id(1)};
    ufline = regexprep(ufline,'.smr.*','');
    ufline = regexprep(ufline,'.smrx.*','');
    ufline = regexprep(ufline,'.*/','');
    ufline = regexprep(ufline,'.*\\','');
    a = regexprep(name,'.*/','');
    a = regexprep(a,'.mat$','');
    if ~strcmp(a,ufline)
        Expt = AddError(Expt, 'Name mismatch %s.smr vs %s.mat',ufline,a);
    end
end
j = j+1;
id = find(strncmp('sb',aText.text,2));
codes(id) = j;
vartypes(j) = 'X';
slens(j) = 0;
id = find(strncmp('#',aText.text,1));
j = j+1;
codes(id) = j;
vartypes(j) = 'X';
slens(j) = 0;

cmid = find(codes == find(strncmp('cm=',findstrs,3)));
rfid = find(strncmp('cm=rf',aText.text(cmid),5));
rfid = cmid(rfid);
bkid = find(strncmp('cm=noback',aText.text(cmid),9));
bkid = cmid(bkid);

if readmethod == 1 && isfield(Trials,'EndTxt')
    endtimes = Trials.EndTxt;
else
    readmethod = 0;
    endtimes = Trials.End;
end
waitloop = round(length(aText.text)/20); %update waitbar 20 times
%waitbar(0,sprintf('parsing %d txt lines',length(aText.text)));
nc = 0;

tendid = aText.codes(:,4) == 2 & aText.codes(:,1) ~= STARTEXPT;% & ismember(aText.codes(:,1),[WURTZOKW WURTZOK BADFIX]);
tendid = tendid';
tdelay = 0;
flipdir = 1;
trialendtime = Trials.End(1);


acodes = aText.codes(:,1);
acodes(find(tendid)) = 100;
newstarts = acodes == FRAMESIGNAL;

%Anything added to  StimTypes here must exist in Stimulsu
StimTypes.char.OptionCode = 1;
Stimulus.OptionCode = '+se';
Stimulus.dfx = 0;

inexpt = 0;



f = {'Seedseq' 'Stimseq' 'Phaseseq' 'cLseq' 'cRseq' 'xoseq' 'yoseq' 'rptframes' 'nsf' 'ntf' 'stimname'};
for j = 1:length(f)
    Stimulus.(f{j}) = {};
    StimTypes.cell.(f{j}) = 1;
end

f = {'st' 'op' 'Fr' 'optionb'};
for j = 1:length(f)
    Stimulus.(f{j}) = 0;
    StimTypes.num.(f{j}) = 1;
end

f = {'Flag' 'explabel' 'exptvars' 'Vn'};
for j = 1:length(f)
    Stimulus.(f{j}) = '';
    StimTypes.char.(f{j}) = 1;
end

nfpj=0;
bsctr = 0;
waitloop = round(length(aText.text)/20); %update waitbar 20 times
E.explabel = '';
electrodestrs = {};
Ignored = [];
state.manual = 0;
%to allow conditional stops to check for processing
%of a given line
xxid = find(strncmp('fp',aText.text,2));

for j = 1:length(aText.text)
%    aText.text{j} =  deblank(aText.text{j});
    txt = aText.text{j};
    t = aText.times(j);
    acode = acodes(j);
    dcode = aText.codes(j,4);
    if codes(j) > 0
        c = codes(j);
        slen = slens(c);
        vartype = vartypes(c);
        ss = txt(1:slen);
        if length(txt) > slen
            if txt(slen+1) == '='
                val= txt(slen+2:end);
            else
                val= txt(slen+1:end);
            end
        else
            val = [];
        end
    else
        id = strfind(txt,'=');
        if ~isempty(id)
            slen = id(1)-1;
            val= txt(slen+2:end);
        elseif length(txt) > 2
            slen = 2;
            val= txt(slen+1:end);
        else
            slen =0;
            val = '';
        end
        if ~isempty(txt)
            ss = txt(1:slen);
        end
        vartype = 'N';
    end
    if slen > 0
        if vartype == 'X'
            if sum(strncmp(txt,{'RightHemi' 'Electrode' ' Electrode' 'Experime'},8)) %lines to include in Comments
                txtid = [txtid j];
                a = binoc.InterpretLine(txt);
                if (isfield(a,'electrode') && ~strcmp(a.electrode,'default')) || ~isfield(Peninfo,'electrode')
                    Peninfo = CopyFields(Peninfo, a);
                    Peninfo.trode = txt;
                end
                if isfield(a,'electrode')
                    electrodestrs{end+1} = a.electrode;
                end
                if isfield(a,'user')
                    Peninfo.user = a.user;
                end
                id = strfind(txt,'Contact');
                if length(id)
                    x = id(1);
                    id = strfind(txt(id:end),' ');
                    Peninfo.probesep = sscanf(txt(id(1)+x:end),'%d');
                end
            elseif strncmp(txt,'Off at',6) %Storage turned off - should be outside trial
                if intrial
                    fprintf('Storage Off in Trial at %.1f',aText.times(j));
                end
            elseif strncmp(txt,'StartDepth',10)
                Stimulus.StartDepth = str2num(txt(11:end));
                StimTypes.num.StartDepth = 1;
            elseif strncmp(txt,'CLOOP',5)
                Stimulus.CorLoop = 1;                
            elseif strncmp(txt,'ch10',4)
%ch10 is R. Everywhere else order is L,R.                
                if txt(5) == '+';
                    Stimulus.eyetest(2) = 1;
                else
                    Stimulus.eyetest(2) = 0;
                end
            elseif strncmp(txt,'ch11',4)
                if txt(5) == '+';
                    Stimulus.eyetest(1) = 1;
                else
                    Stimulus.eyetest(1) = 0;
                end
            elseif strncmp(txt,'rw',2)
                nrw = nrw+1;
                [Trials.rws(nrw), ok] = sscanf(val,'%f');
                Trials.rwset(nrw) = t;
            elseif strncmp(txt,'st',2)
                Stimulus.st = find(strcmp(val, stimnames));
                Stimulus.st = Stimulus.st -1;
                Stimulus.stimname = val;
            elseif strncmp(txt,'mtop=op',7)
                Stimulus.OptionCode = txt(8:end);
                state = CheckOptions(state, Stimulus.OptionCode);
            elseif strncmp(txt,'fl+',3)
                Stimulus.Flag = txt(3:end);
            elseif strncmp(txt,'mtrP=',5)
                Stimulus.Phaseseq = sscanf(txt(6:end),'%d');
                StimTypes.cell.Phaseseq = 1;
                if length(Stimulus.Phaseseq) > length(Stimulus.Stimseq)
                    Stimulus.Phaseseq = Stimulus.Phaseseq(1:length(Stimulus.Stimseq));
                end
                if trial > length(Trials.Start)
                    Trials.Phaseseq{trial} = Stimulus.Phaseseq;
                elseif instim ~= 1  && trial > 1 && t < Trials.Start(trial)
                    Trials.Phaseseq{trial-1} = Stimulus.Phaseseq;
                end
                if Stimulus.id == 6136 || Stimulus.id > 570
                    trial;
                end
            elseif strncmp(txt,'mtFl=',5)
                if length(txt) > 5
                    Stimulus.mtFl = sscanf(txt(6:end),'%d');
                    StimTypes.cell.mtFl = 1;
                end
            elseif strncmp(txt,'mtFi',5) && length(txt) > 50 %in earler versions, this is more reliable
%Becuase char buffer length was running out for mtFn
                Stimulus.mtFi = sscanf(txt(6:end),'%d');
                StimTypes.cell.mtFi = 1;
                skipframes = CheckFrameDiffs(Stimulus.mtFi);
                if ~isempty(skipframes)
                        Trials.rptframes{trial-1} = id;
                        Trials.ndrop(trial-1) = length(id);
                end
            elseif strncmp(txt,'mtFn=',5) && length(txt) > 50
                framets = sscanf(txt(6:end),'%f');
                id = find(diff(framets(1:end-1)) > 1.5);
                if diff(framets(end-1:end)) > 2.8
                    id = cat(1,id ,length(framets));
                end
                if instim ~= 1 && trial > 1
                    if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                        Trials.rptframes{trial-1} = id;
                        Trials.ndrop(trial-1) = length(id);
                        Trials.framet{trial-1} = framets;
                    elseif instim == 2 && t > Trials.End(trial)
                        Trials.framet{trial} = framets;
                        Trials.ndrop(trial) = length(id);
                        Trials.rptframes{trial} = id;
                    end
                    
                end
            elseif strncmp(txt,'mtrX=',5)
                Stimulus.xoseq = sscanf(txt(6:end),'%d');
                if instim ~= 1 && trial > 1
                    if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                        Trials.xoseq{trial-1} = Stimulus.xoseq;
                    elseif instim == 2 && t > Trials.End(trial)
                        Trials.xoseq{trial} = Stimulus.xoseq;
                    end
                    
                end
            elseif strncmp(txt,'mtrY=',5)
                Stimulus.yoseq = sscanf(txt(6:end),'%d');
                if instim ~= 1 && trial > 1
                    if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                        Trials.yoseq{trial-1} = Stimulus.yoseq;
                    elseif instim == 2 && t > Trials.End(trial)
                        Trials.yoseq{trial} = Stimulus.yoseq;
                    end
                end
            elseif strncmp(txt,'mtco=',5)
                Stimulus.Stimseq = sscanf(txt(6:end),'%d');
            elseif strncmp(txt,'mtcL=',5)
                Stimulus.cLseq = sscanf(txt(6:end),'%x');
                StimTypea.cLseq = 'cell';
                if instim ~= 1 && trial > 1
                    if trial > length(Trials.Start) || t < Trials.Start(trial)
                        Trials.cLseq{trial-1} = Stimulus.cLseq;
                    elseif instim == 2 && t > Trials.End(trial)
                        Trials.cLseq{trial} = Stimulus.cLseq;
                    end
                end
            elseif strncmp(txt,'mtcR=',5)
                Stimulus.cRseq = sscanf(txt(6:end),'%x');
                StimTypea.cRseq = 'cell';
                if sum(Stimulus.cRseq < 0) > 1
                    Stimulus.cRseq = sscanf(txt(6:end),'%x');
                end
                if instim ~= 1 && trial > 1
                    if trial > length(Trials.Start) || t < Trials.Start(trial)
                        Trials.cRseq{trial-1} = Stimulus.cRseq;
                    elseif instim == 2 && t > Trials.End(trial)
                        Trials.cRseq{trial} = Stimulus.cRseq;
                    end
                end
            elseif strncmp(txt,'mtrS=',5)
                if aText.codes(j,1) == ENDSTIM
                    istim = 2;
                end
                Stimulus.Stimseq = sscanf(txt(6:end),'%d');
                if Stimulus.id >= 5430
                    Stimulus.Stimseq;
                end
                if instim ~= 1 && trial > 1
                    if trial > length(Trials.Start) || t < Trials.Start(trial)
                        Trials.Stimseq{trial-1} = Stimulus.Stimseq;
                    elseif instim == 2 && t > Trials.End(trial)
                        Trials.Stimseq{trial} = Stimulus.Stimseq;
                    end
                else
                    instim;
                end
            elseif strncmp(txt,'mtse=',5)
                Stimulus.Seedseq = sscanf(txt(6:end),'%d');
                if instim ~= 1 && trial > 1
                    if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                        Trials.Seedseq{trial-1} = Stimulus.Seedseq;
                    elseif instim == 2 && t > Trials.End(trial)
                        Trials.Seedseq{trial} = Stimulus.Seedseq;
                    end
                else
                    instim;
                end
                
            elseif strncmp(txt,'Nf',2) %comes after end TRIAL, not every stim
                if ~instim & trial > 1
                    Trials.Nf(trial-1) = str2num(val);
                end
            elseif strncmp(txt,'mtet=',5)
                id = strfind(txt,'Fr');
                Stimulus.Fr = sscanf(txt(id+3:end),'%d');
            elseif strncmp(txt,'mtxo=',5)
                if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                    lastix = ix;
                    ix = find([Expts.midtrial] > trial);
                    if ix
                        ix = ix(1);
                        Expts(ix).xovals = sscanf(txt(6:end),'%f');
                    end
                end
            elseif strncmp(txt,'dx:',3)
                Trials.dxseq{trial} = sscanf(txt(4:end),'%f');
            elseif strncmp(txt,'mtei=',5)
                if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                    lastix = ix;
                    ix = find([Expts.midtrial] > trial);
                    if ix
                        ix = ix(1);
                        Expts(ix).e1vals = sscanf(txt(6:end),'%f');
                    else
                        err = sprintf('No Expt for mtei at trial %d',trial);
                        Expt = AddError(Expt, err);
                        ix = lastix;
                    end
                end
            elseif strncmp(txt,'mte3=',5)
                if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                    lastix = ix;
                    ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
                    if ix
                        Expts(ix).e3vals = sscanf(txt(6:end),'%f');
                    else
                        ix = lastix;
                    end
                end
            elseif strncmp(txt,'mte2=',5)
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                Expts(ix).e2vals = sscanf(txt(6:end),'%f');
            elseif strncmp(txt,'Off at',5)
            elseif strncmp(txt,'EndExpt',5)
                Stimulus.explabel = '';
                Stimulus.exptvars = '';
                if inexpt
                    ix = ix+1;
                end
                inexpt = 0;
            elseif strncmp(txt,'sonull',5)
            elseif strncmp(txt,'NewConnect',7)
            elseif strncmp(txt,'BGCS Version',7)
            elseif strncmp(txt,'testflag',7)
            elseif strncmp(txt,'rptframes ',3)
                Stimulus.rptframes = sscanf(txt(10:end),'%d');
            elseif strncmp(txt,'seof',4)
                Stimulus.seof = sscanf(txt(5:end),'%d');
            elseif strncmp(txt,'/local',6)
                Stimulus.imprefix = txt;
            elseif strncmp(txt,'imve ',3)
                [a,b] = sscanf(txt,'imve %f,%f %f');
                if length(a) > 2 %otherwise txt is corrupt
                    Stimulus.imver = a(1);
                    Stimulus.imseed = a(2);
                    if length(a) > 2 & a(3) < 1
                        Stimulus.impx = a(3);
                    end
                end
            elseif strncmp(txt,'Sa:',3)
                if aText.codes(j,4) == 1
                    a =  aText.codes(j,1);
                end
            elseif strncmp(txt,'op',2)
                a = sscanf(txt(3:end),'%f,%f');
                if txt(3) == '='  %char description
                    Stimulus.OptionCode = txt(4:end);
                    state = CheckOptions(state, Stimulus.OptionCode);
                elseif ~isempty(a)
                    Stimulus.op = a(1);
                end
                if length(a) > 1
                    Stimulus.optionb = a(2);
                end
                if isempty(Stimulus.op)
                    fprintf('Missing op stim %d\n',trial);
                    Stimulus.op = 0;
                end
            elseif strncmp(txt,'annTyp',6)
                Stimulus.annTyp = sscanf(txt(7:end),'%f');
            elseif ss(1) == '#'
                comment = ss;

            elseif strncmp(txt,'fp',2)
                a = sscanf(val,'%f');
%NB stimulus.dfx records new absolute positions of fp, not changes 
                if length(a) > 1
                    if instim == 1
                        nfpj = nfpj+1;
                        Stimulus.dfx(nfpj) = a(1);
                        Stimulus.dfy(nfpj) = a(2);
                    else
                        Stimulus.fx = a(1);
                        Stimulus.fy = a(2);
                            Stimulus.dfx = a(1);
                            Stimulus.dfy = a(2);
                        nfpj = 0;
                    end
                end
            elseif strncmp(txt,'fx',2) || strncmp(txt,'fy',2)
                %if really in a stimulus, make note of new fx but keep original also
                %if instim ==2, don't change anything.
                a = sscanf(val,'%f');
                if instim == 1
                    Stimulus.(['d' ss]) = a;
                elseif instim == 0
                    Stimulus.(ss) = a;
                    Stimulus.(['d' ss]) = a;
                    lastfix.(ss) = a;
                elseif instim == 2 %Post stim, but trial counter not yet incremented. Store value
                    Stimulus.(ss) = a;
                    Stimulus.(['d' ss]) = a;
                    lastfix.(ss) = a;
                end
            elseif strncmp(txt,'cm=rf',5)
                a = sscanf(txt,'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
                Stimulus.rf = a;
            elseif strncmp(txt,'expname',5)
                if strcmp(val,'NotSet')
                    Stimulus.explabel = '';
                else
                    Stimulus.explabel = val;
                end
            elseif strncmp(txt,'VisualArea',8)
                Stimulus.Vn = val;
            elseif strncmp(txt,'useexvals',8)
                x = sscanf(val,'%d');
                useexvals(1:length(x)) = x;
            else
                Ignored(end+1) = j;
            end %end of 'X'
% Here is where we process standared lines with numeric valcues            
%
%
       elseif aText.codes(j,4) == 2  % this was FROM spike 3
            if acode == 3 && instim == 1 %end stim
                instim = 2;
            elseif acode == STARTEXPT
                inexpt = 1;
                if isfield(Stimulus,'exptvars')
                    Stimulus.exptvars = '';
                end
            end
        elseif acode > 0
            if acode == 5 %stim start
                instim = 1;
                Stimulus.Seedseq = {};  %% these must be set for each stim
                Stimulus.Stimseq = {};
                Stimulus.xoseq = {};
                Stimulus.yoseq = {};
                Stimulus.Phaseseq = [];
                Stimulus.cLseq = [];
                Stimulus.cRseq = [];
                gotend = 0;
                nfpj=0;
                bsctr = bsctr+1;
            elseif acode == 3 %end stim
                instim = 2;
            elseif acode == STARTEXPT 
                inexpt = 1;
                if readmethod == 1
                    id = find([Expts.end] < aText.times(j));
                    ix = length(id)+1;
                end
            elseif acode == ENDEXPT %end stim
                Stimulus.explabel = '';
                if inexpt
                    ix = ix+1;
                end
                inexpt = 0;
            end
        elseif strncmp(txt,'vs',2) || strncmp(txt,'sq',2)
            try
                [a, ok] = sscanf(val,'%f');
                if ~ok
                    Stimulus.(ss) = val;
                else
                    Stimulus.(ss) = a(1);
                    if length(a) > 2
                        Stimulus.FlipDir = a(3);
                    end
                end
            catch
                ok = 0;
                Stimulus.FlipDir = 0;
            end
        elseif strncmp(txt,'EndStim',7) %finished reading all text related to last stim
            gotend = 1;
        elseif strncmp(txt,'manexpt=',8)
            Stimulus.manexpt = txt(9:end);
        elseif sum(strncmp(txt,{'manexvals' 'manexpvals'},8))
            if strncmp(txt,'manexvals',8)
                Stimulus.stimid = sscanf(txt,'manexvals%d');
            else
                Stimulus.stimid = sscanf(txt,'manexpvals%d');
            end
            id = strfind(txt,' ');
            if ~isempty(id)
                a = sscanf(txt(id(1)+1:end),'%f');
            end
            Stimulus.exvals = a;
        elseif strncmp(txt,'exvals',6)
            if txt(7) == ' ' %exvals +space = Manual stim, written by binoc
                a = sscanf(txt,'exvals %f %f %f %d');
                a = sscanf(txt(8:end),'%f');
                if length(a) > 3
                    Stimulus.stimid = a(end);
                end
                Stimulus.exvals = a;
            elseif txt(7) == '='
                a = sscanf(txt(8:end),'%f');
            else
                a = sscanf(txt(7:end),'%f');
            end
            if isfield(Stimulus,'et') && useexvals(1)
                Stimulus.(Stimulus.et) = a(1);
            end
            if length(a) < 2
                if state.manual == 0
                    fprintf('Only 1 expt val in %s\n',txt);
                end
            elseif isfield(Stimulus,'e2') 
                if sum(strcmp(Stimulus.e2,{'backMov' 'Dc'}))
                    Stimulus.(Stimulus.e2) = a(2);
                end
            end
            if length(a) > 2
                Stimulus.ex3val = a(3);
                if isfield(Stimulus,'e3')
                    if strcmp(Stimulus.e3,'mixac')
                        Stimulus.(Stimulus.e3) = a(3);
                    end
                end
            end
            %        elseif strncmp(txt,'id',2)
            %           Stimulus.id = sscanf(txt(3:end),'%d')
        elseif strncmp(txt,'bt',2)
            g = sscanf(txt,'bt%d spkgain %f');
            if length(g) > 1 & g(2) > 1
                Stimulus.SpikeGain = g(2);
            end
            if length(g) > 0
                ExptStartTime = g(1);
            end
        elseif ~isstrprop(ss(1),'alphanum')
            if ~isstrprop(ss(1),'print')
                fprintf('Non-Printing Name %s\n',txt);
            else
                fprintf('Line %d. Not a valid Var Name:%s\n',j,txt);
            end
            
        elseif vartype == 'C'
            if isspace(val(1)) 
                val = val(2:end);
            end
            Stimulus.(ss) = val;
            if ~isfield(StimTypes.char,ss)
                StimTypes.char.(ss) = 'char';
            end
            if strncmp(ss,'et',2)
                Stimulus.(ss) = val;
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                if(ix) Expts(ix).et = val; end
            end
            if strncmp(ss,'e2',2) & length(ix)==1
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                Expts(ix).e2 = val;
            end
            if strncmp(ss,'e3',2) & length(ix) ==1
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                Expts(ix).e3 = val;
            end
             if strncmp(ss,'VisualArea',8)
                Stimulus.Vn = val;
             end

        else
%            if strmatch(ss,{'0' '1' '2' '3' '4'})
            if regexp(ss,'^[0-9]')
                ss = ['x' ss];
            end

            try
            [Stimulus.(ss), ok] = sscanf(val,'%f');
            if ~ok 
                if ~isfield(StimTypes.char,ss)
                    if strncmp('xx',val,2)
                        StimTypes.ignore.(ss) = 1;
                        Stimulus = rmfield(Stimulus,ss);
                    elseif isfield(Stimulus,ss) && isfield(StimTypes.num,ss)%already numeric
                        fprintf('Cannot Make Stimulus field %s a char (%s) - Already Numeric\n',ss,val);                        
                    else
%this is dangerous with arbitray string inputs, where first 2 chars may match an
%old code. e.g. "centerxy"  makes ce a char variable. 
                        StimTypes.char.(ss) = 1;
                        fprintf('Adding %s as char (%s)\n',ss,val);
                        Stimulus.(ss) = val;
                    end            
                else
                    Stimulus.(ss) = val;
                end
            else
                if ~isfield(StimTypes.num,ss) && ~isfield(StimTypes.cell,ss)
                    StimTypes.num.(ss) = 1;
                end
            end
            catch
                ok = 0;
            end
        end
    end
    if acode == FRAMESIGNAL
        newstart = 1;
    else
        newstart = 0;
    end
    if acode == STARTEXPT
        inexpt = 1;
    end
    if tendid(j)  && newstart == 0 %end of data for current trial
        correctdir = 0;
        if isfield(Stimulus,'OptionCode') && ...
                (~isempty(strfind(Stimulus.OptionCode,'+2a')) || ~isempty(strfind(Stimulus.OptionCode,'+afc'))) ...
                & isfield(Stimulus,'vs') || bitand(Stimulus.op,PSYCHBIT);
            [a,b] = max(abs([Stimulus.vs(1) Stimulus.sq(1)]));
%
% historically negative respdir means +ve sacccade value
% for exactly oblique saccades, sign of vertical component does it. 
            if b == 1
                correctdir = -sign(Stimulus.vs(1));
            else
                correctdir = -sign(Stimulus.sq(1));
            end
            Stimulus.rwdir = correctdir;
%            Stimulus.FlipDir = 1;
        end
%        fprintf('t%.1f end %.1f c%.0f trials %d\n',t,endtimes(trial),aText.codes(j,1),trial)
        
%Real ON/Off times are set from the events above. But need to know that
%text following WURTZOK applies to the next stimulus. So instim = 2 means
%text has been received ending trial, but not officially over yet. 
%trial gets incremented at end stim. So response applies to trial -1
%but check that the time is sensible. Can get one of these events when
%storage is off, resetting the last stored trials
       if trial > 1
           tdelay = aText.times(j) - Trials.End(trial-1);
       else
           tdelay  = 0;
       end
       if isfield(Trials,'FlipDir')
           flipdir = Trials.FlipDir(trial-1);
       else
           flipdir = 1;
       end
        if aText.codes(j,1) == WURTZOKW & trial > 1 & tdelay < 10000
            Trials.RespDir(trial-1) = -1 * correctdir.*flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
  %          Trials.score(trial-1) = 0;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == WURTZOK & trial > 1  & tdelay < 10000
            Trials.RespDir(trial-1) = 1 * correctdir.* flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
 %           Trials.score(trial-1) = 1;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == BADFIX
            if trial <= length(Trials.End)
            if aText.times(j) > Trials.End(trial) && settrials== 0
                Trials.End(trial) = aText.times(j);
            end
            if t < Trials.Start(trial) 
%                if trial > 1 && t < Trials.End(trial-1)
%                    Trials.Result(trial-1) = 0;
%                end
            else
                Trials.Result(trial) = 0;
            end
            end
            instim = 2;
        elseif tdelay > 10000 && ismember(aText.codes(j,1),[WURTZOK WURTZOKW])
            if aText.codes(j,3) > 0 
                fprintf('Long delay (%.2f ms) at %.2f between StimON down and es in text\n',tdelay/10,t./10000);
            elseif state.online
                fprintf('Long delay (%.2f ms) at %.2f between StimON down and es in text\n',tdelay/10,t./10000);
            end
        else %shouldb't happen
            trial = trial;
       end
    end
    instims(j) = instim;
    if length(endtimes) >= trial && t >= endtimes(trial) & instim & trial <= ntrials && newstart == 0

%need to read past the end of the last trial a litle way to get things like
%Stimseq which come afterwards
        endlinetimes(trial) = t;
        
        if trial > length(Trials.Start) || trial > length(Trials.End)
            break;
        end
        if length(endtimes) > trial && t > endtimes(trial+1) 
            if sum(endtimes(trial:end)) > 0
            fprintf('?Skipped Trial(%d) at %.1f,%.1f\n',trial,endtimes(trial),t);
            end
        end
        took = (now-tstart) * 24 * 60 *60;
        if took - lasttook > 30
            fprintf('%.0fsec..',took);
            lasttook = took;
        end
        if isfield(Stimulus,'st')
            Stimulus.inexpt = inexpt;
        Stimulus.endevent = aText.codes(j,1);

        Trials = SetTrial(Stimulus, StimTypes, Trials, trial, ntrials);
%reset these parts of Stim in case they changed in Trial
        Stimulus.CorLoop = 0;
        Stimulus.fx = lastfix.fx;
        Stimulus.fy= lastfix.fy;
        Stimulus.rptframes = [];
        if ~isempty(ustimmarkch)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1 ...
                & bitand(ustimmarkch.codes(:,1),1));
            if length(marks)
            Trials.uStimt{trial} = ustimmarkch.times(marks);
            elseif isfield(Trials,'optionb') && bitand(Trials.optionb(trial),64)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1);
            end
            if bitand(Trials.optionb(trial),64)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1);
            end
        end
%        AllStimuli(trial) = Stimulus; %fails when extra element is added to Stimulus
        endtxtid(trial) = j;
        trial = trial+1;
        end
        instim = 0;
    elseif trial > length(Trials.Start)
        if gotend
          break;
        end
% if we have passed stim Start time, text must refer to next stimulus. 
%But don't increment based on codes coming back from Spike
    elseif t > Trials.Start(trial) && instim == 0 && aText.codes(j,4) ~= 2
        instim = 1;
    elseif length(endtimes) >= trial && t >= endtimes(trial)
        if instim
            if newstart
            end
        end
    end
    if guiok && mod(j,waitloop) == 0
        waitbar(j/length(aText.text));
    end
end

if guiok
    delete(waitbar(1));
    drawnow;
end

if state.profiling == 1
    profile viewer;
end

if isfield(Trials,'ve') && iscellstr(Trials.ve)
    Expt.version = Trials.ve{end};
    for j = 1:length(Trials.ve)
        if strncmp(Trials.ve{j},'binoclean',8)
            x = sscanf(Trials.ve{j}(11:end),'%f');
            ve(j) = 10+x(1);
            if length(x) > 1
                ve(j) = ve(j) + x(2)./100;
            end
        else
            ve(j) = 0;
        end
    end
    Trials.ve = ve;
    if length(Trials.ve) < length(Trials.Start)
        Trials.ve(length(Trials.Start)) = median(ve);
    end
end

if length(unique(electrodestrs)) > 1
    Peninfo.electrodestrs = unique(electrodestrs);
end

for j = 1:length(Expts)
    if Expts(j).lasttrial <= Expts(j).firsttrial
        Expts(j).result = 0;
        cprintf('red','Expt %d has no trials\n',j);
    end
end

if  length(Expts) > 1 && isempty(Expts(end).et)
    Expts(end).et = Expts(end-1).et;
    Expts(end).e2 = Expts(end-1).e2;
    Expts(end).e3 = Expts(end-1).e3;
end

if ~isfield(Expts,'result') || isempty(Expts(end).result)
    Expts(end).result = 0;
end

cmid = setdiff(cmid,[rfid bkid]);
cmid = union(cmid,txtid);
Expt.Comments.text = {aText.text{cmid}};
Expt.Comments.times = aText.times(cmid);
Expt.Comments.Peninfo = Peninfo;
Expt.Comments.Peninfo.trode = BuildProbeString(Peninfo);
Expt.Header.CreationDate = CreationDate(aText);

if length(Trials.op) < length(Trials.End) 
        Trials = SetTrial(Stimulus, StimTypes, Trials, length(Trials.Start),ntrials);
end 
if isfield(Trials,'RespDir') &  length(Trials.RespDir) < length(Trials.End)%fill in final trial
        Trials = SetTrial(Stimulus, StimTypes, Trials, length(Trials.Start),ntrials);
end
idx = find(Trials.Result < 0);
if ~isempty(idx)
    fprintf('%d Trials missing End\n',length(idx));
end
fprintf('Text->stims %.2f\n',toc);
tic;



fn = fieldnames(Trials);
ntrials = length(Trials.Start);
cellids = {};
for j = 1:length(fn)
    if iscellstr(Trials.(fn{j}))
        cellids = {cellids{:} fn{j}};
        if length(Trials.(fn{j})) < ntrials
            Trials.(fn{j}){ntrials} = '';
        end
    end
end
eid = GetExptNumber(name);
if eid > 0
    Expt.exptno = eid;
end

[Trials, H] = AddNetFileData(name, Trials, state, Expt.Header);
for j = 1:length(H.errors)
    if H.errortype(j) < 100
        Expt = AddError(Expt,'-noshow',H.errors{j});
    end
end

if isfield(Trials,'mtFi') && H.fixedrpt == 0 %from .bnc file
    for j = 1:length(Trials.Start)
        if j < length(Trials.mtFi) && ~isempty(Trials.mtFi{j})
            Trials.rptframes{j} = CheckFrameDiffs(Trials.mtFi{j});
        end
    end    
end
for j = 1:length(Trials.Start)
    for k = 1:length(cellids)
        if isempty(Trials.(cellids{k}){j})
            Trials.(cellids{k}){j} = '';
        end
    end
end
fprintf('Clearing %d cells took %.2f\n',length(cellids),toc);
trial = 1;
instim = 0;
colors = mycolors;
lastspk = 1;
maxspk = 1;
tpause = 0;
if playspikes
    GetFigure('SpikeV');
end
postdur = 500;

tic;
for trial = 1:length(Trials.Start)
    if trial > length(Trials.End) || isnan(Trials.End(trial)) 
        Trials.End(trial) = Trials.Start(trial) + nomdur;
        Trials = SetTrial(Stimulus, StimTypes, Trials, trial, length(Trials.Start));
    end
    spkids = find(Spks.times > Trials.Start(trial)-preperiod & Spks.times < Trials.End(trial)+postperiod);
    Trials.Spikes{trial} = (Spks.times(spkids) - Trials.Start(trial));
    if spkids
        spkid(trial,:) = [spkids(1) spkids(end)];
        Trials.Cluster{trial} = Spks.codes(spkids,1);
        frames = [];
 %this must be just for caliration of timing. Can't see need for this
 %if we recording real data
        if testframe && ~isempty(frametimes)
            for k = spkids'
                id = find(frametimes < Spks.times(k));
                if id
                    frames = [frames frametimes(id(end))-Spks.times(k)];
                end
            end
        end
        if isempty(frames)
            Trials.Frames(trial) = NaN;
        else
            Trials.Frames(trial) = frames(1);
        end
    else
        spkid(trial,:) = [0 0];
        Trials.Frames(trial) = NaN;
        Trials.Cluster{trial} = [];
    end
    if playspikes
        for spk = spkids; 
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            drawnow;
            hold on;
            energy(spk) = sum(diff(adc).^2);
            svar(spk) = var(adc);
            vw(spk) = svar(spk)/energy(spk);
            if tpause
                pause(tpause);
            end
        end
        %        GetFigure('SpikeXY');
        subplot(2,1,1);
        plot(energy(lastspk:maxspk),vw(lastspk:maxspk),'.');
        hold on;
        %       GetFigure('SpikeV');
        subplot(2,1,2);
        hold off;
    end
end
fprintf('Spikes  take %.3f\n',toc);
if exist(dstimch,'var') & exist('stimchange','var')
    for j = 1:length(stimchange.times)
        et = stimchange.times(j) * 10000;
        id = find(Trials.Start < et);
        if length(id)
            t = id(end);
            if Trials.End(t) > et
                id = find(frametimes < et);
                if ~isempty(id)
                    diffs(j) = et-frametimes(id(end));
                    et = frametimes(id(end));
                end
                Trials.Events{t}{1,2} = et - Trials.Start(t);
                Trials.Events{t}{1,1} = 'ns';
            end
        end
    end
else
    Trials.Events{trial} = [];
end


tic;
trial = 1;
for spk = 1:size(Spks.times)
    t = Spks.times(spk);
    if t > (Trials.End(trial)+postdur) & instim
        Trials.Spikes{trial} = (Spks.times(lastspk:maxspk) - Trials.Start(trial));
        spkid(trial,:) = [lastspk, maxspk];%
%        Trials(trial).Frames = frames;
        Trials.Cluster{trial} = Spks.codes(lastspk:maxspk,1);
        trial = trial+1;
        instim = 0;
        lastspk = spk;
        elseif t > Trials.Start(trial)
            instim = 1;
        end
        if instim & playspikes
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            energy(spk) = sum(diff(adc).^2);
            svar(spk) = var(adc);
            vw(spk) = svar(spk)/energy(spk);
            hold on;
        end
        maxspk = spk;
        if trial > length(Trials.Start)
            break;
        end
end
fprintf('Spikes Took %.3f\n',toc);
fixfile = strrep(idxfile,'idx','fix'); 
if exist(fixfile)
    load(fixfile);
    f = fields(fixTrials);
    [tid, fid] = ismember(fixTrials.id,Trials.id);
    tid = find(tid);
    for j = 1:length(f)
        if ~strcmp(f{j},'id')
            if ~isfield(Trials,f{j})
                Trials.(f{j}) = ones(size(Trials.Start));
            end
            Trials.(f{j})(fid(tid)) = fixTrials.(f{j})(tid);
        end
    end
end
%Trials = Trials(2:end);
AllData.datenum = clusterdate;
AllData.quickload = quickload;
AllData.Events = Events;
AllData.Text = aText;
AllData.Spikes = Spks;
if exist('Ch32','var') && strcmp(Ch32.title,'DigMark')
    Expt.DigMark = Ch32; %So its saved
end

bncname = name;

stored = find(bitand([Trials.op],STOREBIT));
% bitand STOREBIT (16) identifies if storage was on;
%Trials = Trials(stored);

%Expts #20 seems to be empty, but storage was on...
if isfield(Trials,'TrueEnd') & length(Trials.TrueEnd) < length(Trials.End)
    Trials.TrueEnd(length(Trials.End)) = NaN;
end
if isfield(probes,'var')
for j = 1:length(probes)
   probes(j).var = strrep(probes(j).var,'A','');
end
end


for j = 1:length(Expts)
    if Expts(j).lasttrial < Expts(j).firsttrial
        Expts(j).result = 0;
        cprintf('red','Expt %d has no trials\n',j);
    end
end
Expt.Trials = Trials;
Expt.Spkid = spkid;
Expt.Probes = probes;
Expt.bstimes= bstimes;
Expt.estimes= estimes;
Expt.ExptList = Expts;
if state.combineexpts && sum([Expts.result] == 2) > 1
    id = find([Expts.result] ==2);
    for j = 1:length(id)
        tid = Expts(id(j)).firsttrial:Expts(id(j)).lasttrial;
        Expts(id(j)).nf = prctile(Trials.nf(tid),50);
    end
    if state.combineexpts ==2 %automatic: check they are the same
        if length(unique({Expts(id).et})) == 1 && ...
           length(unique({Expts(id).e2})) == 1 && ...
           length(unique([Expts(id).nf])) == 1 && ...
           length(unique({Expts(id).e3})) == 1
              go = 1;
        else
            go = 0;
        end
    else
        go = 1;
    end
    if go
        Expts(id(1)).end = Expts(id(end)).end;
        Expts(id(1)).endi = Expts(id(end)).endi;
        Expts(id(1)).lasttrial = Expts(id(end)).lasttrial;
        [Expts(id(2:end)).result] = deal(20);
    end
end
if exist('mainsch','var') && isfield(mainsch,'times')
    Expt.mainstimes = mainsch.times * 10000;
end
if exist('framech','var') && isfield(framech,'times')
    Expt.frametimes = framech.times * 10000;
end

Header.Name = BuildName(name);
Header.Spike2Version = version;
Header.unstored = nonstore;
Header.CreationDate = CreationDate(aText);
Header.ReadMethod = readmethod;
if isfield(Expt,'DataType')
    Header.DataType = Expt.DataType;
end
clear Text;

if isempty(Expts) %nothing in this file
    fclose(logfid);
    return;
end

if mkidx == 1
        fprintf('Saving index %s\n',idxfile);
        if logfid > 0
        fprintf(logfid, '%s Saving index %s\n',datestr(now),idxfile);
        end
        Expt.bstimes = bstimes;
        Expt.estimes = estimes;
        Expt.Header = Header;
        Expt.state = rmfields(Expt.state,'toplevel');
    save(idxfile,'Expt','Expts');
    iExpts = Expts;
end

%Looked like a Typo.  ".Stimon"  is used when reading index from scratch
%".SetStim" is used for SetStimulusVals
%Also wrong arguments. Looks like this was just an earlier attempt at
%fixing that was not deleted
Trials = AddStimValsFromFile(strrep(name,'.mat','.SetStim'),Trials);

if length(Trials.Result) > length(Trials.Start)
    Trials.Result = Trials.Result(1:length(Trials.Start));
end
if isfield(Expts,'firsttrial')
    if timeoffset
        args = {args{:} 'timeoffset' timeoffset};
    end
[Expts, Expt] = SortExpts(Expts, Trials, Header, thecluster, Expt, state);
Expts = AddComments(Expts,Expt);
if isempty(Expts)
    nc = sum([Expt.ExptList.result]==19);
    cprintf('blue','No Expts. %d were cancelled. %d not saved\n',nc,unsavedexpts);
end
end

if mkidx == 1
    tExpts = Expts;
    ExptList = MkExList(Expts);
    Expts = iExpts;
    Expt.starttrial = starttrial;
    save(idxfile,'Expt','Expts','ExptList');
    Expts = tExpts;
end
if saveexpts
    SaveExpts(name, Expts);
end
if mkidx
%    save(idxfile,'Expt','Expts');
end
fclose(logfid);
WriteErrors(idxfile, Expt);


function  StimCh = AddStimsFromFile(AddTxtFile, StimCh)
%Add Stim On/Off events in pairs: format
%timeon timeoff
%timeon -1  DELETES a bogus event (quick on/off pairs with no matching
%serial events
%
%File name mnkXXX.Stimon

fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%f %f','delimiter','\n'); %ON, Off
    for j = 1:length(a{2})
        id = find(StimCh.times < a{1}(j));
        if length(id) < length(StimCh.times) 
            id = id(end);
            if a{2}(j) > 0
                lvl = cat(1,StimCh.level(1:id), 1, 0 ,StimCh.level(1+id:end));
                t = cat(1,StimCh.times(1:id), a{1}(j), a{2}(j), StimCh.times(1+id:end));
                StimCh.level = lvl;
                StimCh.times = t;
            else %delete pair of events after a(1)
                StimCh.level(id+[1:2]) = [];
                StimCh.times(id+[1:2]) = [];
            end
        else %just add end stim
            StimCh.level(end+1) = 0;
            StimCh.times(end+1) = a{2}(j);
        end
    end
    fclose(fid);
end

function  Events = AddEventsFromFile(AddTxtFile, Events)
%Add Events. Each line of the file should contain 
%time code
%time is in second, not ticks
%File name mnkXXX.AddEvents

codes = zeros(1,size(Events.codes,2));
fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%f %f','delimiter','\n'); %ON, Off
    for j = 1:length(a{1})
        id = find(Events.times < a{1}(j));
        if isempty(id)
            id = 0;
        end
        if length(id) < length(Events.times)
            id = id(end);
            t = cat(1,Events.times(1:id), a{1}(j), Events.times(1+id:end));
            if a{2}(j) == -1001;
                codes(1) = 255;
                codes(3) = 1;
            else
                codes(1) = a{2}(j);
            end
            allcodes = cat(1,Events.codes(1:id,:), codes ,Events.codes(1+id:end,:));
            Events.codes = allcodes;
            Events.times = t;
            fprintf('Adding Event %d at %.3f\n',a{2}(j),a{1}(j));
        end
    end
    fclose(fid);
end

function Trials = AddStimValsFromFile(name, Trials)
fid = fopen(name,'r');
if fid > 0
    mycprintf('blue','Reading Stimulus Properties from %s\n',name);
    a = textscan(fid,'id%d %s','delimiter','\n');
    tid = a{1};
    s = a{2};
    for j = 1:length(a{1})
        if strncmp('badexpt',s{j},7) %error in file - ignore this expt
            Trials.baddexpts(tid(j)) = 1;
        else
        id = find(Trials.id  == a{1}(j));
        if length(id) ==1
            x = sscanf(s{j}(4:end),'%f');
            if strncmp(s{j},'dx:',3)
                Trials.dxvals{id} = x;
            elseif strncmp(s{j},'ce:',3)
                Trials.cevals{id} = x;
            end
        end
        end
    end
end

function [Expt, Expts] = LoadExptsFile(name)
    if exist(name)
        suffix = 0;
        load(name);
        Expt.state.nospikes = 1;
        Expt = CopyFields(Expt,Tidx);
        id = regexp(name,'\.[0-9]*Expts.mat');
        if ~isempty(id)
            suffix = sscanf(name(1+id(1):end),'%d');
        end
        for j = 1:length(Expts)
            Expts{j}.Header.loadname = strrep(name,'Expts.mat','.mat');
            Expts{j}.Header.suffix = suffix;
        end
        
    else
        Expt = [];
        Expts= {};
    end

function SaveExpts(name, Expts)
%make separate files for each expt on disk so that can access in parallel
%when needed
for j = 1:length(Expts)
    outfile = strrep(name,'.mat',['Expt' num2str(j) '.mat']);
    Expt = Expts{j};
    save(outfile,'Expt');
end



function [T, H] = AddNetFileData(name, T, state, Header)

recalc = state.buildbnc;
H.params = {};
H.errors = {};
H.fixedrpt =0;
BinocFile = regexprep(name,'\.[0-9]*.mat','.bnc');
BinocFile = regexprep(BinocFile,'([0-9]*).mat','$1.bnc');
BinocBFile = regexprep(name,'\.[0-9]*.mat','.bno'); %manually created version from online
BinocBFile = regexprep(BinocBFile,'([0-9]*).mat','$1.bno'); %manually created version from online
if ~exist(BinocFile)
    BinocFile = BinocBFile;
    
end
BinocMatFile = regexprep(name,'\.[0-9]*.mat','.bnc.mat');
BinocMatFile = regexprep(BinocMatFile,'([0-9]+).mat','$1.bnc.mat');
if isempty(strfind(BinocMatFile,'.bnc')) %failed to make name - probably online. Don't overwrite!
    return;
end

if state.online
    if state.needbnc == 0
        return;
    end
    [a,b,c,d] = GetMonkeyName(BinocFile);
    BinocFile = strrep(BinocFile,'\','/');
    x = regexprep(BinocFile,'online/.*','');
    trueBinocFile = sprintf('%s/data/%s/%s/%s%s.bnc',x,a,c,a,c);
    if exist(trueBinocFile)
        BinocFile = trueBinocFile;
    end
end

%of BinocExpt is stored in the app data, then its date should be the same
%ad the .mat file (BinocFile), so just test that

if exist(BinocFile) %check for changes in .bnc after making .mat
    d = dir(BinocFile);
    a = dir(BinocMatFile);
    if length(a) ~= 1 || length(d) ~= 1  || d.datenum > a.datenum
        fprintf('Re-Reading %s\n',BinocFile);
        recalc = 1;
    end
end
loaded = 0;
%Allow BinocExpt to be stored as appdata in 0, in case calling man times without GUI
%BinocExpt.readdate is the date when the .bnc file was read
%.loaddate records when the .mat file was loaded
if isfield(state,'toplevel') && (isfigure(state.toplevel) || double(state.toplevel) ==0) && isappdata(state.toplevel,'BinocExpt') && recalc == 0
    BinocExpt = getappdata(state.toplevel,'BinocExpt');
    fprintf('Using cached BinocExpt (%s)\n',datestr(BinocExpt.readdate));
    loaded = 1;
elseif isappdata(0,'BinocExpt') && recalc == 0
    B = getappdata(0,'BinocExpt');
    if ~isfield(B,'loadname')
        fprintf('Cached BinocExpt No loadname');
        B.loadname = B.filename;
    end
    if strcmp(B.loadname,BinocMatFile) && exist(BinocMatFile) && B.loaddate > a.datenum
        loaded = 1;
        BinocExpt = B;
        loaded = 1;
        fprintf('Using BinocExpt from Root Window (%s)\n',datestr(BinocExpt.readdate));
    end
end
if loaded == 0
if exist(BinocFile) && (~exist(BinocMatFile) || recalc == 1)
    
    fprintf('Reading %s\n',BinocFile);
    BinocExpt = ReadSerialTrials(BinocFile);
    fprintf('Saving %s...',BinocMatFile);
    ts = now;
    save(BinocMatFile,'BinocExpt');
    fprintf('Took %.2f\n',mytoc(ts));
    BinocExpt.loadname = BinocMatFile;
    BinocExpt.loaddate = now;
elseif exist(BinocMatFile)
    ts = now;
    fprintf('Reading %s',BinocMatFile);
    load(BinocMatFile);
    fprintf('Took %.2f\n',mytoc(ts));
    BinocFile =  BinocMatFile;
    BinocExpt.loadname = BinocFile;
    BinocExpt.loaddate = now;
else 
    return;
end
end
if ~isempty(strfind(BinocBFile,'.bno')) && exist(BinocBFile) && ~isfield(BinocExpt,'bnofile') %botha .bnc and a .bno.  Used to fix missing things like mtFl    
    BinocBExpt = ReadSerialTrials(BinocBFile);
    f = setdiff(fields(BinocBExpt.Trials),fields(BinocExpt.Trials));
    for j = 1:length(BinocExpt.Trials)
        tid = find([BinocBExpt.Trials.id]== BinocExpt.Trials(j).id);
        Tr = BinocBExpt.Trials(tid);
        if isempty(tid)
            for k = 1:length(f)
                BinocExpt.Trials(j).(f{k}) = Tr.(f{k});
            end
        end
        
    end
    newid = setdiff([BinocBExpt.Trials.id],[BinocExpt.Trials.id]);
    newt = find(ismember([BinocBExpt.Trials.id],newid));
    for j = 1:length(newt)
        BinocExpt.Trials = CopySFields(BinocExpt.Trials,0,BinocBExpt.Trials(newt(j)));
    end
    if ~isempty(BinocBExpt.cellfields)
        f = fields(BinocBExpt.cellfields);
        for j = 1:length(f)
            BinocExpt.cellfields.(f{j}) = 2;
        end
    end
    if ~isempty(BinocBExpt.vectorfields)
        f = fields(BinocBExpt.vectorfields);
        for j = 1:length(f)
            BinocExpt.vectorfields.(f{j}) = 2;
        end
    end
    BinocExpt.bnofile = BinocBFile;
    if ~isfield(BinocExpt,'loadname')
        BinocExpt.loadname = BinocMatFile;
    end
    BinocExpt.bnofields = f;
    save(BinocMatFile,'BinocExpt');
    if isfield(state,'toplevel') && isfigure(state.toplevel)
        setappdata(state.toplevel,'BinocExpt',BinocExpt);
    end

end

name = GetFileName(BinocFile);
T.netfields = {};
sf = fields(T);
nf = fields(BinocExpt.Trials);
netf = nf;
f = setdiff(nf, sf);
smrid = T.id;
bT = BinocExpt.Trials;
netid = [bT.id];
[id, sid, nid]  = intersect(smrid, netid);

if isfield(state,'toplevel') && isfigure(state.toplevel)
    setappdata(state.toplevel,'BinocExpt',BinocExpt);
else
    setappdata(0,'BinocExpt',BinocExpt);
end

if ~isempty(BinocExpt.cellfields)
cf = fields(BinocExpt.cellfields);
H.params = {H.params{:} cf};
nbt = [];
%If there are dupicate ids in .bnc, find closest match
%based on calendar times.  N.B. this relies on clock syncrhonization
%between machines.  Assuming any error is smaller than time separation
%bewteen repeats of id, this should be OK. 
resolved = ones(1,length(nid));
for j = 1:length(nid)
    nbt(j) = sum(netid == netid(nid(j)));
    if nbt(j) > 1 && isfield(BinocExpt.Trials,'TrialTime') %duplicate Ids in bnc file       
        t  = expt.ConvertTime(Header,T.Start(sid(j)));
        tid = find(netid == netid(nid(j)));
        [a,b] = min(abs(t - [BinocExpt.Trials(tid).TrialTime]));
        if a < 1  %can't be wrong day!!
            nid(j) = tid(b);
        else
            resolved(j) = 0;
        end
    end
end

if sum(nbt>1)
    fprintf('%d Duplicate Ids in %s, %d unresolved\n',sum(nbt>1),BinocFile, sum(resolved ==0));
end

cf = intersect(cf,fields(bT));
for j = 1:length(cf)
    nx = [];
    for k = 1:length(nid)
        nx(k) = length(bT(nid(k)).(cf{j}));
    end
    if max(nx) > 1 %only convert ->cell array if there are mulitple values for this expt
        if isfield(T,cf{j}) && ~iscell(T.(cf{j}))
            x = T.(cf{j});
            T.(cf{j}) = {};
            for k = 1:length(x) %first convert existing field to cell
                T.(cf{j}){k} = x(k);
            end
        end
        for k = 1:length(nid)
            
            nf = T.nf(sid(k))+1;
            if isfield(bT,'Nf') && ~isempty(bT(nid(k)).Nf) && bT(nid(k)).Nf > 0
                
                nf = bT(nid(k)).Nf; %actuall number painted
            end
            if isfield(T,'Fn') & T.Fn(sid(k)) > 1
                Fn = T.Fn(sid(k));
                nf = ceil(nf/Fn) .* Fn;
            end
            x = bT(nid(k)).(cf{j});
            xf = 0;
            if length(x) > nf && T.Result(sid(k)) > 0
                if length(x) > nf+xf
                    H.errors{end+1} = sprintf('%d Frames in %s, but only %d shown id%d',length(x),name,nf,id(k));
                    x = x(1:nf);
                    H.errortype(length(H.errors)) = 1;
                    if T.inexpt(sid(k)) ==0
                        H.errortype(length(H.errors)) = 101;
                    end
                        
                end
            end
            T.(cf{j}){sid(k)} = x;
        end
        if max(sid) < length(T.Start)
            T.(cf{j}){length(T.Start)} = [];
        end
    else
        for k = 1:length(nid)           
            nf = T.nf(sid(k))+1;
            x = bT(nid(k)).(cf{j});
            if ~isempty(x)
                if ~isfield(T,cf{j})
                elseif iscell(T.(cf{j}))
                    T.(cf{j}){sid(k)} = x;
                else
                    T.(cf{j})(sid(k)) = x;
                end
            end
        end
        
    end
end
T.netfields = {T.netfields{:} cf{:}};

%now check for dropped frames
for k = 1:length(nid)
        xf = 0;
        if isfield(bT,'mtFn')
            if ~isempty(strfind(T.OptionCode{sid(k)},'+FN'))
                H.fixedrpt = H.fixedrpt+1;
                T.rptframes{sid(k)} = bT(nid(k)).mtFl;
                skips = bT(nid(k)).mtFl;
                xf = length(skips);
                
%nskip should report how many fraems the skipped frame was up for

                [skipframes, errval, nskip] = CheckFrameDiffs(diff(bT(nid(k)).mtFn),'frameperiod',1);
                skipid = skipframes < length(bT(nid(k)).mtFn)-1;
                T.skipframes{sid(k)} = skipframes(skipid);
                T.nskipped{sid(k)} = nskip;
                skipid = find(skipid);
                if xf == 1 && ~isempty(skipframes) && skipframes(1) - skips(1) > 100
                    T.rptframes{sid(k)} = []; %no real skips
                    xf = 0;
                elseif xf > 0 && ~isempty(skipframes)
                    fprintf('Frame Skip %d,%d    ',skips(1),skipframes(1));
                    rptf = [];
                    for f = 1:length(skipid)
                        ns = nskip(skipid(f));
                        rptf(end+1:end+ns) = skipframes(skipid(f));
                    end
                    T.rptframes{sid(k)} = rptf;
                elseif ~isempty(T.skipframes{sid(k)})
                    fprintf('Frame Skip at %d (not mtFl)',T.skipframes{sid(k)}(1));
                    T.rptframes{sid(k)} = skipframes;
                end
            else
                T.rptframes{sid(k)} = bT(nid(k)).mtFn;
                xf = length(bT(nid(k)).mtFn);
            end
        end
        if xf> 0
            fprintf('Dropped %d frames id%d\n',xf,id(k));
        end
    end
end

checkf = intersect({'skipframes' 'nskipped'},fields(T));
for j = 1:length(checkf)
    if length(T.(checkf{j})) < length(T.Start)
        T.(checkf{j}){length(T.Start)} = [];
    end
end
%fields that should come from bnc
checkf = intersect({'expvars'},netf);
for j = 1:length(checkf)
    for k = 1:length(nid)
        T.(checkf{j}){sid(k)} = bT(nid(k)).(checkf{j});
    end
end

if ~isempty(BinocExpt.vectorfields)
cf = fields(BinocExpt.vectorfields);
H.params = {H.params{:} cf};
for j = 1:length(cf)
    if ~isfield(T,cf{j})
        for t = 1:length(nid)
            T.(cf{j}){sid(t)} = bT(nid(t)).(cf{j});
        end
    elseif ~iscell(T.(cf{j}))
    end
end
T.netfields = {T.netfields{:} cf{:}};
end

if ~isempty(H.errors)
    [H.errors,id] = unique(H.errors);
    H.errortype = H.errortype(id);
end
fprintf('Finished Reading %s\n',BinocFile);


function state = CheckOptions(state, str)

if strfind(str,'+exm')
    state.manual = 1;
else
    state.manual = 0;
end

function str = BuildProbeString(P)

str = '';
if isfield(P,'electrode')
    str = [str 'Electrode ' P.electrode ' '];
end
if isfield(P,'tube')
    str = [str 'Tube ' P.tube ' '];
end
if isfield(P,'hemishpere')
    str = [str P.hemisphere 'HemiSphere '];
end

function ix = FindExptn(Expts, readmethod, bsctr, trial, ix)
if readmethod == 1
    trial = bsctr;
end
    
if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
    lastix = ix;
    ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
    if ix
        ix = ix(end);
    else
        ix = lastix;
    end
end

function ExptList = MkExList(Expts)
ExptList = [];
if ~iscell(Expts)
    return;
end
for j =1:length(Expts)
        ExptList(j).expname = Expts{j}.Header.expname;
        ExptList(j).start = Expts{j}.Header.Start;
        ExptList(j).end = Expts{j}.Header.End;
        ExptList(j).et = Expts{j}.Stimvals.et;
        ExptList(j).e2 = Expts{j}.Stimvals.e2;
        ExptList(j).e3 = Expts{j}.Stimvals.e3;
end

function WriteErrors(idxfile, Idx)

if isfield(Idx,'errs') & length(Idx.errs)
    ename = strrep(idxfile,'idx.mat','err.txt');
    fid = fopen(ename,'a');
    fprintf(fid,'%s\n',Idx.errs{:});
    fclose(fid);
end
if isfield(Idx,'errdata') && ~isempty(Idx.errdata)
    errfile = BuildFileName(idxfile,'error');
    if ~exist(errfile)
        X.Errors = Idx.errdata;
    else
        X = load(errfile);
        if iscell(X.Errors)
            X.Errors = CellToStruct(X.Errors);
        end
        ne = length(X.Errors);
        f = fields(Idx.errdata);
        for j = 1:length(Idx.errdata)
            for k = 1:length(f)
                X.Errors(ne+j).(f{k}) = Idx.errdata(j).(f{k});
            end
        end
    end
    save(errfile,'-struct','X')
end


function [aText, Text, newlines] = AddText(AddTxtFile, aText, Text)
%Add text reads lines of the form
%t txt
%where t is in seconds, not timestamps
%Addts these to the txt record.
%
%t = -1000 puts lines at the end
%t  = -1 finds all lines that begin with the same code 
%        and replaces them
%can also add Stimulus Values by id with
%idxxxx str
%t delete
%          removes line # t
%t replace xxxx puts all text following 'replace' into text element t
%
% N.B  Add.txt is read separately above, looking for StarExpt/EndExpt
% additions


newlines = 0;
[a,b,c] = fileparts(AddTxtFile);
if ~strcmp(c,'.txt')
    return;
end
SpkDefs;
fid = fopen(AddTxtFile,'r');
ts = now;

if fid > 0
    cprintf('blue','Adding Text from %s\n',AddTxtFile)
    a = textscan(fid,'%d %s','delimiter','\n');
    if isempty(a{1})
        if isempty(aText.text)
            return;
        end
        b = textscan(fid,'%2s%d %s','delimiter','\n');
        frewind(fid);
        idid = find(strncmp('id',aText.text,2));
        for j = 1:length(idid)
            tid(j) = sscanf(aText.text{idid(j)},'id%d');
        end
        ti = b{2};
        for j = 1:length(ti)
            id = find(tid == ti(j));
            if ~isempty(id)
                t(j) = aText.times(idid(id(end)));
            else
                t(j) = NaN;
            end
        end
        s = b{3};        
    else
    t = a{1};
    s = a{2};
    end
    fclose(fid);
    
    newlines = 0;
    did = find(strncmp(s,'delete',6));
    if ~isempty(did) && ~isempty(aText.times)
        for j = 1:length(did)
            cprintf('blue','%d deleting %s\n',t(j),aText.text{t(did(j))});
        end
        id = setdiff(1:length(aText.times), t(did));
        aText.text = aText.text(id);
        aText.times = aText.times(id);
        aText.codes = aText.codes(id,:);
    end
    for j = 1:length(s)
        %?why do we look for whitespace?? removed Aug 2010.
        if(strncmp(s{j},'replace',6))
            id = findstr(s{j},' ');
            aText.text{t(j)} = s{j}(id(1)+1:end);
        else
            id = findstr(s{j},' ');
            id = [];
            if length(id)
                txt = s{j}(id(1)+1:end);
            else
                txt = s{j};
            end
            if j < 500
                cprintf('blue','%d Adding Text %s\n',t(j),s{j})
            end
            if ~isempty(aText.times)
                if t(j) < 0  && t(j) > -1000%special case for fixing lines
                    if strncmp(s{j},'cm=rf',5)
                        id = find(strncmp('cm=rf',aText.text,5));
                        for k = 1:length(id)
                            aText.text(id(k)) = s(j);
                        end
                    end
                elseif t(j) == -1000  %put these at the end
                    aText.text = {aText.text{:} txt};
                    aText.times = [aText.times; max(aText.times)+1];
                    aText.codes = [aText.codes; 0 0 0 0];
                elseif t(j) == 0 || t(j) <= aText.times(1)
                    aText.text = {txt aText.text{:}};
                    aText.times = [0; aText.times];
                    aText.codes = [0 0 0 0; aText.codes];
                elseif sum(strcmp(s(j),{'delete' 'badexpt'})) %special lines not going into text
                    txt = '';
                else
                    aText.text{end+1} = txt;
                    aText.times(end+1)=t(j);
                    aText.codes(end+1,:)=0;
                    if strcmp(txt,'EndExpt')
                        aText.codes(end,1) = ENDEXPT;
                    end
                    newlines = newlines+1;
                end
            end
            if strncmp(s{j},'VisualArea',8) && t(j) == 0
                id = find(strncmp('VisualArea',aText.text,8));
                for k = 1:length(id)
                    aText.text{id(k)} = s{j};
                end
            end
            if ~isempty(txt)
                Text.text(end+1,1:length(txt)) = txt;
                newlines = newlines+1;
            end
        end
    end
    if newlines && ~isempty(aText)
        [t, tid] = sort(aText.times);
        aText.text = aText.text(tid);
        aText.times = aText.times(tid);
        aText.codes = aText.codes(tid,:);
    end
    fprintf('Took %.2f\n',mytoc(ts));
end



function Expts = AddComments(Expts, Expt)

if isfield(Expt,'Comments')

for j = 1:length(Expt.Comments.times)
    id = find(Expt.Trials.Start < Expt.Comments.times(j));
    if isempty(id)
        Expt.Comments.id(j) = Expt.Trials.id(1);
    else
        Expt.Comments.id(j) = Expt.Trials.id(id(1));
    end
end

ts = 0;
    for j = 1:length(Expts)
        te = Expts{j}.Header.trange(2)+100000;
        if j > 1
            ts = min([ ts Expts{j-1}.Header.trange(2)]);
        end
        if j < length(Expts)
            te = max([ te Expts{j+1}.Header.trange(1)]);
        end
        cid = find(Expt.Comments.times > ts & Expt.Comments.times < te);
        if ~isempty(cid)
            bid = find(strncmp('cm=back=',Expt.Comments.text(cid),8));
            for k = 1:length(bid)
                %these lines are before start expt. If after end its for
                %next expt
                if Expt.Comments.times(cid(bid(k))) < Expts{j}.Header.trange(2)
                    Expts{j} = ParseExptComment(Expts{j}, Expt.Comments.text{cid(bid(k))});
                end
            end
        end
        Expts{j}.Comments.text = {Expt.Comments.text{cid}};
        Expts{j}.Comments.times = [Expt.Comments.times(cid)];
        cid = find(strncmp('cm=VisualArea',Expt.Comments.text,12));
        for k = 1:length(cid)
            Expts{j}.Header.Area{k} = Expt.Comments.text{cid(k)}(15:end);
        end
        if isfield(Expt.Comments,'Peninfo')
            Expts{j}.Comments.Peninfo = Expt.Comments.Peninfo;
        end
        ts = te;
    end
end



function Trials = SetTrial(Stimulus, StimTypes, Trials, trial, ntrials)


if trial > length(Trials.Start)
    return;
end
fn = fieldnames(Stimulus);
if isfield(Stimulus,'Ro') && isfield(Stimulus,'dx')
    ca = cos(Stimulus.Ro * pi/180);
    sa = sin(Stimulus.Ro * pi/180);
% need to recheck sign conventions here in replay...    
    Stimulus.Op = Stimulus.yo .* ca - Stimulus.xo .* sa;
    Stimulus.Pp = Stimulus.yo .* sa + Stimulus.xo .* ca;
    Stimulus.dO = Stimulus.dx .* sa - Stimulus.dy .* ca;
    Stimulus.dP = Stimulus.dy .* sa + Stimulus.dx .* ca;
end

fnum = fieldnames(StimTypes.num);
fchar = fieldnames(StimTypes.char);
fcell = fieldnames(StimTypes.cell);

fnon = setdiff(fn, [fnum' fchar' fcell']);
trialf = fieldnames(Trials);

%for fields in Stimulus, not in trials, and not cell or char
% pr allocate memory
newf = setdiff(fn,[trialf' fchar' fcell']);
for k = 1:length(newf)
    Trials.(newf{k})(1:ntrials,1) = NaN;
end    

for k = 1:length(fchar)
    F = fchar{k};
    if ~isfield(Trials,F) || ischar(Trials.(F)) || iscell(Trials.(F))
        if ~isempty(Stimulus.(F))
            Trials.(F){trial} = Stimulus.(F);
        end
    elseif ~isempty(Stimulus.(F)) %if its empty, the test above fiale
        fprintf('%s is Char in Stimulus, not in Trials\n',F);
    end
end

for k = 1:length(fnum)
    Trials.(fnum{k})(trial,1:length(Stimulus.(fnum{k}))) = Stimulus.(fnum{k});
end
for k = 1:length(fnon)
    Trials.(fnon{k})(trial,1:length(Stimulus.(fnon{k}))) = Stimulus.(fnon{k});
end
for k = 1:length(fcell)
    Trials.(fcell{k}){trial} = Stimulus.(fcell{k});
end

if isempty(Stimulus.st) || isnan(Stimulus.st) || isnan(Trials.st(trial))
    fprintf('Unrecognized stimulus %s trials %d\n',Stimulus.stimname,trial);
    Stimulus.st
end
if isfield(Trials,'rwset') && isfield(Trials,'rws')
id = find(Trials.rwset < Trials.Start(trial));
eid = find(Trials.rwset < Trials.End(trial));
if ~isempty(id)
    Trials.rw(trial) = Trials.rws(id(end));
end
end
if isfield(Trials,'RespDir') && length(Trials.RespDir) < length(Trials.Start)
    Trials.RespDir(length(Trials.Start)) = 0;
%    Trials.FlipDir(length(Trials.Start)) = 1;
end


function Header = SetHeader(Header, AllTrials, igood, f)
%Set Header field using value of a cell string in selected trials
%Allows for empty Cells, which unique does not

    if isfield(AllTrials,f)
        if iscell(AllTrials.(f))
%if its a cellstr, but some are emppy, this messes up counts.
%so restrict this to non-empty elements
            C = AllTrials.(f)(igood(igood < length(AllTrials.(f))));
            [ok, b] = CellToMat(C);
            [a, b] = Counts(C(ok > 0));
            if ~isempty(b) && ~isempty(a)
                Header.(f) = b{1};
            elseif isfield(Header,f)
                Header = rmfield(Header,f);
            end
        else
            Header.(f) = unique(AllTrials.(f)(igood));
        end
    end

function [Expts, Idx, state] = SortExpts(AllExpts, AllTrials, Header, thecluster, Idx,state,  varargin)
SpkDefs;
timeoffset = 0;
%stimnames = {'None', 'Gabor', 'RDS', 'Grating', 'bar', 'circle', 'rectangle', 'test', 'square', 'probe', '2grating', 'Cylinder', 'twobar', 'rls', 'annulus', 'rdssine', 'nsines'};
Expts = [];
state.tt = TimeMark(state.tt,'Start Sorting');
spikid = Idx.Spkid;
if ~isfield(Idx,'newerrs')
    Idx.newerrs = 0;
end
if isfield(Header,'frameperiod')
    frameperiod = Header.frameperiod;
    if isnan(Header.frameperiod) && state.needframes == 0
        if isfield(AllTrials,'fz')
            frameperiod = 10000./(median(AllTrials.fz));
        else
            frameperiod = 167;
        end
        Header.frameperiod = frameperiod;
    end
else
    frameperiod = 167;
end

if isfield(Header,'Name') && ~isempty(regexp(Header.Name,'[0-9]\.[0-9]*\.mat'))
    Header.bysuffix = 1;
else
    Header.bysuffix = 0;
end


id = find(diff(AllTrials.id) < 0);
for j = 1:length(id)
    offset = AllTrials.id(id(j))-AllTrials.id(id(j)+1);
    AllTrials.id(id(j)+1:end) = AllTrials.id(id(j)+1:end)+offset;
    AllTrials.idoffset(id(j)+1:length(AllTrials.id)) = offset;
    Idx = AddError(Idx, 'Adding Id Offset %d at ID %d',offset,id(j)+1);
end


findtrial = 0;
if state.alltrials
    usebadtrials = 1;
else
usebadtrials = 0;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'alltrials',5)
        usebadtrials = 1;
    elseif strncmpi(varargin{j},'findtrial',5)
        j = j+1;
        findtrial = varargin{j};
    elseif strncmpi(varargin{j},'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
    end
    j = j+1;
end

outname = strrep(Header.loadname,'.mat','Expts.mat');
if exist(outname) && state.resort == 0
    X = load(outname,'ExptState');
    if isfield(X,'ExptState') && X.ExptState.alltrials == state.alltrials %otherwise need to reload
        fprintf('Using Expts saved in %s\n',outname);
        [Expt, Expts] = LoadExptsFile(outname);
        return;
    end
end

if isfield(AllTrials,'ve') && length(AllTrials.ve) < length(AllTrials.Start)
        AllTrials.ve(length(AllTrials.Start)) = median(AllTrials.ve);
end


for j = 1:length(AllExpts)
    if ~isfield(AllExpts,'result') || isempty(AllExpts(j).result)
        AllExpts(j).result = -1;
    end
end
nalltrials = length(AllTrials.Start);

AllExpts = AllExpts([AllExpts.result] >= 0);
STIM_BAR = 4;
STIM_GRATING = 3;
for j = 1:length(stimnames)
    if find(strcmp(stimnames{j},'bar'))
        STIM_BAR = j-1; %starts with 0
    end
end

fn = fieldnames(AllTrials);
% make a list of fileds that are NOT automatically set 
state.tt = TimeMark(state.tt,'Fixed Names');
ids = find(strcmp('Spikes',fn));
ids = [ids find(strcmp('OptionCode',fn))];
ids = [ids find(strcmp('Start',fn))];
ids = [ids find(strcmp('End',fn))];
ids = [ids find(strcmp('EndTxt',fn))];
%ids = [ids find(strcmp('Result',fn))]; %keep result now can mark badfix
%after stim
ids = [ids find(strcmp('op',fn))]; % has op and optionb
ids = [ids find(strcmp('Stimseq',fn))];
ids = [ids find(strcmp('Seedseq',fn))];
%ids = [ids find(strcmp('rwtimes',fn))];
ids = [ids find(strcmp('xoseq',fn))];
ids = [ids find(strcmp('yoseq',fn))];
ids = [ids find(strcmp('rptframes',fn))];
ids = [ids find(strcmp('framet',fn))];
ids = [ids find(strcmp('cLseq',fn))];
ids = [ids find(strcmp('cRseq',fn))];
ids = [ids find(strcmp('Phaseseq',fn))];
ids = [ids find(strcmp('Cluster',fn))];
ids = [ids find(strcmp('Frames',fn))];
ids = [ids find(strcmp('StartEv',fn))];
ids = [ids find(strcmp('Events',fn))];
ids = [ids find(strcmp('rws',fn))];
ids = [ids find(strcmp('endelay',fn))];
ids = [ids find(strcmp('rwset',fn))];
ids = [ids find(strcmp('estimes',fn))];
ids = [ids find(strcmp('uStimt',fn))];
ids = [ids find(strcmp('Flag',fn))];
ids = [ids find(strcmp('bsstimes',fn))];
ids = [ids find(strcmp('esstimes',fn))];
ids = [ids find(strcmp('ex3val',fn))];
ids = [ids find(strcmp('bsdelay',fn))];
ids = [ids find(strcmp('stimname',fn))];
ids = [ids find(strcmp('explabel',fn))];
ids = [ids find(strcmp('exptvars',fn))];
ids = [ids find(strcmp('imprefix',fn))];
ids = [ids find(strcmp('netfields',fn))];
ids = [ids find(strcmp('pe',fn))];
ids = [ids find(strcmp('psychfile',fn))];
ids = [ids find(strcmp('VisualArea',fn))];
ids = [ids find(strcmp('mo',fn))];

ids = [ids find(strcmp('Trw',fn))];
ids = [ids find(strcmp('uf',fn))];
%ids = [ids strmatch('imver',fn)];
%ids = [ids strmatch('imse',fn)];
state.tt = TimeMark(state.tt,'Set Fields');


%do not include PhaseSeq here - has special cases
seqstrs = {'dxvals' 'cevals'};
cf = {'rwtimes' 'nsf' 'ntf' 'skipframes' 'Vn'};
    cellfields = {};
for j = 1:length(cf)
    f = cf{j};
    if isfield(AllTrials,f) && iscell(AllTrials.(f))  %need to make this general
        cellfields = {cellfields{:} f};
    end
end
seqvars = {};
for j = 1:length(seqstrs)
    id = find(strcmp(seqstrs{j},fn));
    if ~isempty(id)
        ids = [ids id];
        seqvars = {seqvars{:} fn{id}};
    end
end

%exptvars is list of variables manipulated by a manual expt
exptvars = {};
    if isfield(AllTrials,'uStimt') && length(AllTrials.uStimt) < length(AllTrials.Start)
                    AllTrials.uStimt{length(AllTrials.Start)} = [];
    end

fn = fn(setdiff([1:length(fn)],ids));
for nf = 1:length(fn)
    if iscell(AllTrials.(fn{nf}))
    else
        if length(AllTrials.(fn{nf})) < length(AllTrials.Start)
            fprintf('Forcing end values for %s\n',fn{nf});
            AllTrials.(fn{nf})(length(AllTrials.Start)) = 0;
        end
    end
end
state.tt = TimeMark(state.tt,'Checked Ends');
%non-zero values of falsestart should indicate the time gap
if isfield(AllTrials,'TrueStart') & ~isfield(AllTrials,'FalseStart')
    AllTrials.FalseStart = AllTrials.TrueStart;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 0)) = 1;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 1)) = 0;
end
if ~isfield(AllTrials,'Nf') || length(AllTrials.Nf) < length(AllTrials.Start)
    AllTrials.Nf(length(AllTrials.Start)) = NaN;
end

if isfield(AllTrials,'Phaseseq')
    for j = 1:length(AllTrials.Phaseseq)
        if ~isempty(AllTrials.Phaseseq{j})
            id = find(AllTrials.Phaseseq{j} == 0);
            AllTrials.Phaseseq{j}  = 1;
        end
    end
end
if ~isfield(AllTrials,'Events') || length(AllTrials.Events) < length(AllTrials.Start)
    AllTrials.Events{length(AllTrials.Start)} = [];
end

if length(AllTrials.Result) < length(AllTrials.Start)    
    Idx = AddError(Idx, 'Ex %d Result only %d/%d Trials\n',GetExptNumber(Header.Name),length(AllTrials.Result),length(AllTrials.Start));
    AllTrials.Result(length(AllTrials.Start)) = 0;
end
%If there are a minority of trials with False Starts (usually missing the
%StimON channel, or a long delay to the matching FRAMESIGNAL on the serial line
%set the result to 2, so that these can be exluded if necessary.
% only do this if the result is already not 0 (otherwise Add in Badfix
% trials
if isfield(AllTrials,'FalseStart') && sum(AllTrials.FalseStart > 0) < length(AllTrials.FalseStart)/5
    if length(AllTrials.FalseStart)> length(AllTrials.Result)
        AllTrials.FalseStart = AllTrials.FalseStart(1:length(AllTrials.Result));
    end
    AllTrials.Result(find(AllTrials.FalseStart > 0 & AllTrials.Result > 0)) = 2;
end
    %fn is now a list of fields that are set for each trial automatically.
nexpts = 1;
phasevals  = [0:360];
phasevals(2:4) = [pi pi/4 3*pi/4];
lphasevals = [0 pi pi 0]; %%see binoc.c SetRandomPhase();
rphasevals = [0 pi 0 pi];
for nx = 1:length(AllExpts)
    Idx.exptno = nexpts;
    state.tt = TimeMark(state.tt,sprintf('Expt %d',nx));
    clear Trials;
    
    frpt = 1;
    a = AllExpts(nx).firsttrial:AllExpts(nx).lasttrial;
    a = a(a <= nalltrials);
    if isfield(AllTrials,'sM') && median(AllTrials.sM(a)) == 26
        badtr = 1;
    else
        badtr = usebadtrials;
    end
    
    if badtr
        Header.usebadtrials = 1;
        igood = 1:length(a);
        fn = {fn{:} 'Result'};
        ngood = sum(AllTrials.Result(a) > 0);
    else
        igood = find(AllTrials.Result(a) > 0);
        ngood = length(igood);
    end
    ntotaltrials = length(a);
    nt = length(igood);
    igood = a(igood);
    nu = 0; %number of ustim pulses
    % even when using all trials, don't include expts that have no good ones
    if ngood <= 3 && AllExpts(nx).result ~= CANCELEXPT && state.showerrs
        err = sprintf('%s Expt %d only %d/%d good trials',Header.Name,nx,ngood,ntotaltrials);
        Idx = AddError(err, Idx, 0);
    end
    if ~isfield(AllExpts,'result')
        AllExpts(1).result = 1;
    end
    if ~isfield(AllExpts,'e3')
        AllExpts(1).e3 = 'e0';
    end
    if ngood> 3 & igood(1) < length(AllTrials.Start) & ismember(AllExpts(nx).result,[2 0])
       spkids = [];
       needfields = {};
       allframes = [];
       if isempty(AllExpts(nx).e3)
           fprintf('Empty Type Expt %d, trials %d - %d\n',nx, AllExpts(nx).firsttrial, AllExpts(nx).lasttrial)
       elseif sum(strcmp(AllExpts(nx).e3,'ar'))
           needfields = {needfields{:} 'wi' 'hi'};
       end
       if isfield(AllTrials,'optionb')
           AllTrials.optionb(isnan(AllTrials.optionb)) = 0;
       end

       for j = 1:nt
            if findtrial & AllTrials.Start(igood(j)) > findtrial
                findtrial = 0;
            end
            [Trials(j).Start] = AllTrials.Start(igood(j));
            [Trials(j).TrialStart] = AllTrials.Start(igood(j));
            if isfield(AllTrials,'TrueEnd') & AllTrials.TrueEnd(igood(j)) > 0
                [Trials(j).End] = AllTrials.TrueEnd(igood(j));
            else
                [Trials(j).End] = AllTrials.End(igood(j));
            end
            Trials(j).dur = Trials(j).End(end)-Trials(j).Start(1);
            if isfield(AllTrials,'optionb') 
                [Trials(j).uStim] = bitand(AllTrials.optionb(igood(j)),64);
            end
            [Trials(j).op] = AllTrials.op(igood(j));
            [Trials(j).Trial] = AllTrials.Trial(igood(j));
            [Trials(j).id] = AllTrials.id(igood(j));
            if ~isempty(AllTrials.Events{igood(j)})
                Trials(j).Events = AllTrials.Events{igood(j)};
            end
            if isfield(Trials,'uStimt')
                nu = nu + legnth(Trials(j).uStimt);
            end
%            spkids = [spkids spkid(igood(j),1):spkid(igood(j),2)];
        end
 % cant do this. for some reason max spkids) is > size(Spikes.times)
%        Expt.gui.spks = unique(spkids);
%  problem is that we don't rebuilds Spkid if set s probe != 1, so save
%  Spkid no good
 %      
        fastseq = 0;
        
        for nf = 1:length(cellfields)
            if iscellstr(AllTrials.(cellfields{nf}))
                nv = unique({AllTrials.(cellfields{nf}){igood}});
            else
                if length(AllTrials.(cellfields{nf})) < max(igood)
                    AllTrials.(cellfields{nf}){max(igood)} = [];
                end
                [x, nv] = Counts(AllTrials.(cellfields{nf})(igood));
            end
            if length(nv) > 1
                needcellfield(nf) = 1;
            else
                needcellfield(nf) = 0;                
            end            
        end
 rcvars = {};
 
        for nf = 1:length(fn)
            fsz= size(AllTrials.(fn{nf}));
            flen = 1;
            needf = 0;
            if iscell(AllTrials.(fn{nf}))
                if strcmp(fn{nf},'uStimt')
                else
                    Expt.Stimvals.(fn{nf}) = AllTrials.(fn{nf}){a(1)};
                    if iscellstr(AllTrials.(fn{nf}))
                        nv = unique({AllTrials.(fn{nf}){igood}});
                    else
                        [x, nv] = Counts(AllTrials.(fn{nf})(igood));
                        if length(x) ==1 && length(nv) > 1 && isnumeric(nv) %one element in each, but different values
                            nv = {nv};
                        elseif length(x) ==1
                            nv = 1;
                        end
                    end
                end
            else
                if fsz(1) > 1 && fsz(2) > 1
%for now only check unique of element 1. If 1s all same, and 2s all same, but 1 ~= 2, 
%don't need it
                    flen = fsz(2);
                    if fsz(1) < max(igood)
                        AllTrials.(fn{nf})(max(igood),1) = 0;
                    end
                    nv = unique([AllTrials.(fn{nf})(igood,1)]);
                    for k = 2:fsz(2)
                        nvb{k-1} = unique([AllTrials.(fn{nf})(igood,k)]);
                        if sum(~isnan(unique([AllTrials.(fn{nf})(igood,k)]))) >1
                            needf = 1;
                        end
                    end
                else
                    nv = unique([AllTrials.(fn{nf})(igood)]);
                    nvb{1}  = NaN;
                end
            end
            if strcmp(fn{nf},'exvals')
                    for j = 1:nt
                        Trials(j).exvals = AllTrials.exvals(igood(j),:);
                    end
            elseif isnumeric(nv) && (sum(~isnan(nv))> 1 || needf)
                if fsz(1) > 1 && fsz(2) > 1
                    for j = 1:nt
                        Trials(j).(fn{nf}) = AllTrials.(fn{nf})(igood(j),:);
                    end
                else
                    for j = 1:nt
                        Trials(j).(fn{nf}) = AllTrials.(fn{nf})(igood(j));
                    end
                end
                if strcmp(fn{nf},'Nf')
                    if sum(nv) == max(nv) % only 1 value + 0;
                        Expt.Stimvals.(fn{nf}) = max(nv);
                    end
                    Trials = rmfield(Trials,'Nf');
                end
                Expt.Stimvals.(fn{nf}) = prctile([AllTrials.(fn{nf})(igood)],50);
            elseif isnumeric(nv)
                if sum(~isnan(nv)) > 0
                    Expt.Stimvals.(fn{nf}) = nv(~isnan(nv));
                    for k = 2:flen
                        bvals = nvb{k-1}(~isnan(nvb{k-1}));
                        if length(bvals) == 1
                            Expt.Stimvals.(fn{nf})(k) = bvals(1);
                        elseif isempty(bvals)
                            fprintf('No values for %s(%s)\n',fn{nf},k);                            
                        else
                            Expt.Stimvals.(fn{nf})(k) = bvals(1);
                            fprintf('>1 value for %s(%d)\n',fn{nf},k);                            
                        end
                    end
                else
                    Expt.Stimvals.(fn{nf}) = NaN;
                end
            elseif sum(strcmp(fn{nf},cellfields)) %don't put these in rcvar
            elseif length(nv) > 1 %cell arary with more than one value
                vals = [];
                for j = 1:nt
                    Trials(j).(fn{nf}) = AllTrials.(fn{nf}){igood(j)};
                    vals = cat(1, vals(:),AllTrials.(fn{nf}){igood(j)}(:));
                end
                Expt.Stimvals.(fn{nf}) = median(vals);
                rcvars = {rcvars{:} fn{nf}};
            elseif iscellstr(nv) %same string on every trial
                Expt.Stimvals.(fn{nf}) = nv{1};
            elseif iscell(nv) %single cell = single element per trial, but different values
                vals = [];
                for j = 1:nt
                    Trials(j).(fn{nf}) = AllTrials.(fn{nf}){igood(j)};
                    vals = cat(1, vals(:),AllTrials.(fn{nf}){igood(j)}(:));
                end
            end
        end
        if isfield(AllTrials,'dfx') && size(AllTrials.dfx,2) > 1
            for j = 1:nt
                id = find(AllTrials.dfx(igood(j),:) ~= 0);
                Trials(j).dfx = AllTrials.dfx(igood(j),id);
                Trials(j).dfy = AllTrials.dfy(igood(j),id);
            end
            
        end
        if isfield(AllTrials,'rf');
            Expt.Stimvals.rf = median(AllTrials.rf(igood,:));
            Expt.Stimvals.st = mode(AllTrials.st(igood));
        end
    duration = mean([Trials.End] - [Trials.Start]);
    if isnumeric(Expt.Stimvals.et)
        Expt.Stimvals.et = NumberToExptName(Expt.Stimvals.et, Expt);
    end
    if isnumeric(Expt.Stimvals.e2)
        Expt.Stimvals.e2 = NumberToExptName(Expt.Stimvals.e2, Expt);
    end
    et = Expt.Stimvals.et;
    e2 = Expt.Stimvals.e2;
    Header = SetHeader(Header, AllTrials, igood, 'explabel');
    Header = SetHeader(Header, AllTrials, igood, 'exptvars');
    Header = SetHeader(Header, AllTrials, igood, 'idoffset');
    if ~isempty(rcvars)
        Header.rcvars = rcvars;
    end
    if isfield(Header,'exptvars') && ~isempty(Header.exptvars)
        exptvars = split(Header.exptvars,',');
        exptvars = setdiff(exptvars,{'NotSet'});
    end
    if ~isfield(Expt.Stimvals,'e3')
        Expt.Stimvals.e3 = 'e0';
    end
    if isfield(Expt.Stimvals,'Fr') && Expt.Stimvals.Fr > 0
        frpt = Expt.Stimvals.Fr;
    end
    if isfield(AllExpts(nx),'e1vals') & ~isempty(AllExpts(nx).e1vals)
        Expt.e1vals = AllExpts(nx).e1vals;
    end
    if isfield(AllExpts(nx),'xovals') & ~isempty(AllExpts(nx).xovals)
        Expt.xovals = AllExpts(nx).xovals;
    elseif isfield(Expt,'xovals')
        Expt = rmfield(Expt,'xovals');
    end
    if sum(strncmp('ce',{et e2},2))
        if iscell(AllTrials.ce)
            Expt.Stimvals.ce = median(abs(cat(1,AllTrials.ce{igood})));
        else
        Expt.Stimvals.ce = median(abs(AllTrials.ce(igood)));
        end
    end
    if ~isfield(Trials,et) && ischar(et)
        needfields = {needfields{:} et};
    end
    missingfields = setdiff(needfields,fields(AllTrials));
    if ~isempty(missingfields)
    end
    needfields = intersect(needfields,fields(AllTrials)); %Can't add field if not in alltrials

     if sum(strcmp(et,{'Op' 'Pp'})) & isfield(Expt.Stimvals,'rf')
            ca = cos(Expt.Stimvals.rf(5) * pi/180);
            sa = sin(Expt.Stimvals.rf(5) * pi/180);
            rOp = Expt.Stimvals.rf(2) .* ca - Expt.Stimvals.rf(1) .* sa;
            rPp = Expt.Stimvals.rf(2) .* sa + Expt.Stimvals.rf(1) .* ca;
            Expt.Stimvals.rOp = rOp;
            Expt.Stimvals.rPp = rPp;
    if isfield(AllExpts(nx),'e1vals')
% don't mess with the values for interleaved extras
            sid = find(AllExpts(nx).e1vals > -1000);
            if strcmp(et,'Op')
                Expt.e1vals(sid) = AllExpts(nx).e1vals(sid) +rOp;
            elseif strcmp(et,'Pp')
                Expt.e1vals(sid) = AllExpts(nx).e1vals(sid) +rPp;
            end
        end
    end
    if isfield(AllExpts(nx),'e2vals') && ~isempty(AllExpts(nx).e2vals)
        Expt.e2vals = AllExpts(nx).e2vals;
        if strcmp(Expt.Stimvals.e2,'ce') && max(Expt.e2vals) > 1
            bid = find(Expt.e2vals > 1);
            err = sprintf('Fixing ce values (%s->1)',sprintf('%.1f ',Expt.e2vals(bid)));
            Idx = AddError(Idx, err);
            Expt.e2vals(bid) = 1;
        end

    end
                if isfield(Expt.Stimvals,'Fs')
        xovals = [-16:16] .* Expt.Stimvals.Fs;
        yovals = [-16:16] .* Expt.Stimvals.Fs;
        end
    allev = [];
    serrid = [];
    nframes = [];
    durs = [];
    psychtrial = 0;
    seqtrial = 0;
    crtrial = 0; %count # with contrast reversal
%Count trials with these options set and set a flag in header if its moat    
    checkoptions = {'+cr' '+x2' '+exm' '+FN'};  
    optioncounts = zeros(1,length(checkoptions));
    timesexpt = 0;
    nu = 0; %count trials with uStimt;
    for k = 1:nt
        Idx.t = AllTrials.Start(igood(k));
        if ~isempty(AllTrials.Cluster{igood(k)})
        sid = find(ismember(AllTrials.Cluster{igood(k)}, thecluster));
        Trials(k).Spikes = round(AllTrials.Spikes{(igood(k))}(sid));
        Trials(k).count = sum(find(Trials(k).Spikes > 500 & Trials(k).Spikes < duration+500));
        else
            Trials(k).Spikes = [];
            Trials(k).count = 0;
        end
        
        for c = 1:length(exptvars)
            if ~isfield(AllTrials,exptvars{c})
                fprintf('Missing Expt var %s\n',exptvars{c});
            elseif iscell(AllTrials.(exptvars{c}))
                Trials(k).(exptvars{c}) = AllTrials.(exptvars{c}){igood(k)};
            else
                Trials(k).(exptvars{c}) = AllTrials.(exptvars{c})(igood(k));
            end
        end
        stimframes = AllTrials.nf(igood(k));
        Trials(k).sz = AllTrials.wi(igood(k));
        Trials(k).OptionCode = AllTrials.OptionCode{igood(k)};
        if isfield(AllTrials,'Flag')
        Trials(k).Flag = AllTrials.Flag{igood(k)};
        if strfind(Trials(k).Flag,'+mm')
            Trials(k).flatsurf = 1;
        else
            Trials(k).flatsurf = 0;
        end
        end
        if isfield(AllTrials,'uStimt')
        Trials(k).uStimt = round(AllTrials.uStimt{(igood(k))});
        nu = nu + length(Trials(k).uStimt);
        end
        for f = 1:length(needfields)
            Trials(k).(needfields{f}) = AllTrials.(needfields{f})(igood(k));
        end
        if ~isempty(strfind(Trials(k).OptionCode,'+2a')) || ~isempty(strfind(Trials(k).OptionCode,'+afc'))
            psychtrial = psychtrial+1;
        elseif bitand(Trials(k).op,PSYCHBIT)
            psychtrial = psychtrial+1;
        end
        if strfind(Trials(k).OptionCode,'+FN')
            allframes(k) = 1;
        else
            allframes(k) = 0;
        end

        if isfield(AllTrials,'Dc')
            Dcval = AllTrials.Dc(igood(k));
        else
            Dcval = 0;
        end
%Trials.(et) can be missing if et is undef, or if et is ''        
%Test for Dcval was making some manual fS expts get caught
        if strfind(Trials(k).OptionCode,'+fS') & Dcval < 1 & isempty(strfind(Trials(k).OptionCode,'+exm'))
            seqtrial = seqtrial+1;
            fastseq = 1;
            manseq = 0;
        elseif isfield(Trials,et) && length(Trials(k).(et)) > Expt.Stimvals.nf/2 %Manual Fast seq
            seqtrial = seqtrial+1;
            manseq = 1;
            fastseq = length(Trials(k).(et));
%If  set manual  sequence longer than nf, then it stops at nf.            
            if fastseq > stimframes +1 && stimframes > 1
                Idx = AddError(Idx, '-noshow', 'Frame Sequence %d, nf only %d at %.2f id%d\n',fastseq,stimframes,Trials(k).Start(1),Trials(k).id);
            end
        else
            fastseq = 0;
            manseq = 0;
        end
        for c = 1:length(checkoptions)
            if strfind(Trials(k).OptionCode,checkoptions{c})
                optioncounts(c) = optioncounts(c)+1;
            end
        end
        if ~isnan(AllTrials.Frames(igood(k)))
            Trials(k).Frames = AllTrials.Frames(igood(k));
        end
        if bitand(LMONOC,Trials(k).op)
            Trials(k).me  = -1;
        elseif bitand(RMONOC,Trials(k).op)
            Trials(k).me  = 1;
        else
            Trials(k).me  = 0;
        end
        if isfield(Trials,'Fr')
            frpt = Trials(k).Fr;
        end
        
        %If just changing im seed slowly, don't mess wiht end/start
        if frpt > 1 && fastseq == 0
            frpt = 1;
        end
% for image seqs, Stimseq records the order of seeds, and there will not
% be a conversion to stimulus type
%NB. Seedseq is NOT in seqvars, and is handled separately. Don't want to
%build substpace maps, to don't make Start/End into vectors
       needseq = 0;
       for ns = 1:length(seqvars)
           if length(AllTrials.(seqvars{ns}){igood(k)}) > 1
               Trials(k).(seqvars{ns}) = AllTrials.(seqvars{ns}){igood(k)};
               needseq = length(Trials(k).(seqvars{ns}));
           end
       end
        if isfield(AllTrials,'Seedseq') && length(AllTrials.Seedseq{igood(k)}) > 1
            Trials(k).Seedseq = AllTrials.Seedseq{igood(k)};
        end
        for nf = 1:length(cellfields) %fileds that are cell arrays
            f = cellfields{nf};
            if needcellfield(nf) && isfield(AllTrials,f) && length(AllTrials.(f)) >= k
                Trials(k).(f) = AllTrials.(f){igood(k)};
            end
        end
        if manseq == 0 && isfield(AllTrials,'Stimseq') && length(AllTrials.Stimseq{igood(k)}) > 1 && ...
                ((isfield(Expt,'e1vals') && ~isempty(Expt.e1vals)) || isfield(Expt,'xovals')) 
            if AllTrials.Stimseq{igood(k)}(end) > length(Expt.e1vals) %last value sometimes junk               
               AllTrials.Stimseq{igood(k)} = AllTrials.Stimseq{igood(k)}(1:end-1);
            end
            if length(Expt.e1vals) && ~isempty(AllTrials.Stimseq{igood(k)}) && ...
                ((max(AllTrials.Stimseq{igood(k)}) < length(Expt.e1vals) && seqtrial > k/2) || strcmp(et,'backMov'))
            evid = AllTrials.Stimseq{igood(k)}+1;
            if frpt > 1
                evid = evid(1:frpt:end);
            end
            Trials(k).st = ones(size(evid)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(evid)) * Expt.Stimvals.ce;
            Trials(k).me = ones(size(evid)) * Trials(k).me(1);
            if sum(strcmp(et,'Dc'))
                ev = Expt.e1vals(evid);
                Trials(k).(et) = AllTrials.Dc(igood(k));
                Trials(k).(Expt.Stimvals.e2)= ev;
                Trials(k).(e2)(find(ev == ISIGNALFRAME)) = AllTrials.(e2)(igood(k)); %% blanks
                if strcmp(Expt.Stimvals.e2,'or')
                    Trials(k).ori = AllTrials.(e2)(igood(k));
                end
            elseif sum(strcmp(et,'backMov'))
                ev = AllTrials.Stimseq{igood(k)}+1;
                Trials(k).(et) = ev;
             %   eb = Expt.e2vals(AllTrials.Stimseq{igood(k)}+1);
              %  Trials(k).(Expt.Stimvals.e2)= eb;
            elseif sum(strcmp(et,{'ic' 'pR'})) && Dcval > 0 && Dcval < 1
                ev = Expt.e1vals(evid);
                Trials(k).ori = AllTrials.or(igood(k));
                Trials(k).or = ev;
            else
                ev = Expt.e1vals(evid);
                Trials(k).(et) = ev;
                 if isfield(Expt,'e2vals') && length(Expt.e2vals) > max(AllTrials.Stimseq{igood(k)})
                    eb = Expt.e2vals(evid);
                    Trials(k).(Expt.Stimvals.e2)= eb;
                 else
                     eb = zeros(size(ev));
                end
            end
            if isfield(AllTrials,'xoseq') & length(AllTrials.xoseq{igood(k)}) > 0
                Trials(k).xo = xovals(AllTrials.xoseq{igood(k)}+1);
            end
            if isfield(AllTrials,'yoseq') & length(AllTrials.yoseq{igood(k)}) > 0
                Trials(k).yo = yovals(AllTrials.yoseq{igood(k)}+1);
            end
            if isfield(AllTrials,'cLseq') & length(AllTrials.cLseq{igood(k)}) > 0
                Trials(k).cL = AllTrials.cLseq{igood(k)};
                Trials(k).cL(Trials(k).cL < 0) = 0;
                Trials(k).cL(Trials(k).cL > 512) = 0;
            end
            if isfield(AllTrials,'cRseq') & length(AllTrials.cRseq{igood(k)}) > 0
                Trials(k).cR = AllTrials.cRseq{igood(k)};
                Trials(k).cR(Trials(k).cR < 0) = 0;
                Trials(k).cR(Trials(k).cR > 512) = 0;
            end
            if isfield(AllTrials,'Phaseseq') & length(AllTrials.Phaseseq{igood(k)}) > 1
                if Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_BAR
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                elseif Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_GRATING
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                else
                Trials(k).ph = phasevals(AllTrials.Phaseseq{igood(k)}+1);
                end
            end
            Trials(k).Start = Trials(k).Start + [0:length(ev)-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            Trials(k).ce(find(ev == IUNCORR)) = 0; %% uncorr
            Trials(k).me(find(ev == ILEFTMONOC)) = -1; 
            Trials(k).me(find(ev == IRIGHTMONOC)) = 1;
                if size(Trials(k).(et),2) > 1
                    size(evid)
                end
                if isfield(Trials,e2) & size(Trials(k).(e2),2) > 1
                    size(evid)
                end
%in the Lopos X Ropos expt, could exceed max stim combinatinos, then mtei
%and mte2 strings only gave the inidividual Lpos/Rpos values, not all
%combinations
        elseif isfield(AllTrials,'Stimseq') && AllTrials.ve(end) < 4.85 && ...
                ~isempty(AllTrials.Stimseq{igood(k)}) && ...
                fastseq && ...
             max(AllTrials.Stimseq{igood(k)}) > length(Expt.e1vals)
         extras = sum(Expt.e1vals < -999);
         evid = AllTrials.Stimseq{igood(k)}-extras;
            Trials(k).st = ones(size(evid)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(evid)) * Expt.Stimvals.ce;
            Trials(k).me = ones(size(evid)) * Trials(k).me(1);
         n2 = length(Expt.e2vals)-extras;
         n1 = length(Expt.e1vals)-extras;
         e1 = mod(evid,n2)+extras+1;
         id = find(evid <0);
         e1(id) = evid(id)+extras+1;
         eb = floor(evid./length(Expt.e1vals))+1+extras;
         eb(id) = evid(id)+extras+1;
         Trials(k).(et) = Expt.e1vals(e1);
         Trials(k).(et)(id) = 0;
         ev = Expt.e1vals(e1);
         if ~isempty(Expt.e2vals)
         Trials(k).(e2) = Expt.e2vals(eb);
         end
            Trials(k).Start = Trials(k).Start + [0:length(ev)-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            if Trials(k).id == 6138
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            end
            if isfield(AllTrials,'Phaseseq') & length(AllTrials.Phaseseq{igood(k)}) > 0
                if Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_BAR
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                elseif Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_GRATING
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                else
                Trials(k).ph = phasevals(AllTrials.Phaseseq{igood(k)}+1);
                end
            end
        elseif isfield(AllTrials,'Stimseq') && length(AllTrials.Stimseq{igood(k)}) > 2 && seqtrial > k/2
            if strcmp(Expt.Stimvals.et,'serange')
                if isfield(AllTrials,'xoseq') & length(AllTrials.xoseq{igood(k)}) > 0
                    Trials(k).xo = xovals(AllTrials.xoseq{igood(k)}+1);
                end
                if isfield(AllTrials,'yoseq') & length(AllTrials.yoseq{igood(k)}) > 0
                    Trials(k).yo = yovals(AllTrials.yoseq{igood(k)}+1);
                end
            end
            serrid = [serrid k];
            end
        elseif needseq
            Trials(k).Start = Trials(k).Start + [0:needseq-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;            
        end
        if isfield(AllTrials,'rptframes') && ~isempty(AllTrials.rptframes{igood(k)})
            Trials(k).rptframes = AllTrials.rptframes{igood(k)};
        end
        
        durs(k) = Trials(k).End(end) - Trials(k).Start(1);
        if fastseq && manseq %count #stim in manual fast seqs
            nframes(k) = fastseq;
        else
            nframes(k) = length(Trials(k).Start);
        end
    end %% End of looping through trials
    if isfield(Trials,'ob') && sum(cat(1,Trials.ob) < 0) && ~isfield(Trials,'or')
        aid = find([Trials.ob] >= 0);
        bid = find([Trials.ob] < 0);
        [Trials(aid).or] = deal(Expt.Stimvals.or);
        if Expt.Stimvals.or > 25
            [Trials(bid).or] = deal(Expt.Stimvals.or-90);
        else
            [Trials(bid).or] = deal(Expt.Stimvals.or+90);
        end
        for j = 1:length(Trials)
            Trials(j).ob = abs(Trials(j).ob);
        end
    end
    if isfield(Trials,'CorLoop')
        id = find([Trials.CorLoop] ~= 0);
        if ~isempty(id)
            Header.excluded.CorLoop = [Trials(id).id];
        end
        id = find([Trials.CorLoop] == 0);
        if length(id)
            Trials = Trials(id);
            nframes = nframes(id);
            durs = durs(id);
        end
    end
    if isfield(Trials,'flatsurf')
        fs = unique([Trials.flatsurf]);
        if length(fs) == 1
            Trials = rmfield(Trials,'flatsurf');
            Expt.Stimvals.flatsurf = fs;
        end
    end

%some expts may have nfr > 2 on some trials, nut not other
%e.g. image tuning exps (seedseq ->nfr > 2) with interleaved blanks.
%only trim for missing sequence if its really an fastseq expt
    if seqtrial > nt/2 && mean(nframes) > 2
        fastseq = 1;
    else
        fastseq = 0;
    end
    if fastseq && state.alltrials == 0
       id = find(nframes > 1);
       bid = find(nframes  ==1);
       if ~isempty(bid)
                  Header.excluded.noRC = [Trials(bid).id];
       end
       Trials = Trials(id);
       if length(bid) > 2
           Idx = AddError(Idx, 'Ex %d Removing %d Trials becuase no RC sequence\n',nexpts,length(bid));
       end
    end
    if nu == 0  && isfield(Trials,'uStimt')
        Trials = rmfield(Trials,'uStimt');
    end
    if isfield(Trials,'uStim') && sum([Trials.uStim]) == 0
            Trials = rmfield(Trials,'uStim');
    end
    Header.idrange = minmax([AllTrials.id(igood)]);
    if isfield(Trials,'inexpt') && sum([Trials.inexpt] ==0 > 0)
       bid = find([Trials.inexpt] == 0);
       Header.excluded.noexpt = [Trials(bid).id];
       id = find([Trials.inexpt] > 0);
       ids = unique([Trials(bid).id]);
       if ~isempty(bid)
           Idx.t = Trials(bid(1)).Start(1);
           if length(ids) <= 1
               Trials = Trials(id);
               Idx = AddError(Idx, 'Ex %d Removing %d Trials becuase not in Expt\n',nexpts,length(bid));
           else
               Idx = AddError(Idx, 'Ex %d Has  %d Trials not in Expt\n,',nexpts,length(bid));
           end
       end
    end

    id = [];
    if isfield(Trials,'delay')
        id = find(~isnan([Trials.delay]));
    end
    if length(id) > length(Trials)/2
        nid = find(isnan([Trials.delay]));
        if isfield(Idx,'frametimes')
            for j = 1:length(nid)
                nft = sum(Idx.frametimes > Trials(nid(j)).Start(1) & Idx.frametimes < Trials(nid(j)).End(end)); 
                if nft == 0
                    cprintf('red','No Frame Signals for Trial Id %d - storge probably off\n',Trials(nid(j)).id);
                end
            end
        end
        if length(nid)
            if isfield(Trials,'FalseStart')
            fprintf('Delay Nan at %.2f(%.2f), id%.0f\n',Trials(nid(1)).Start(1),Trials(nid(1)).FalseStart, Trials(nid(1)).id);
            else
            fprintf('Delay Nan at %.2f, id%.0f\n',Trials(nid(1)).Start, Trials(nid(1)).id);
            end
        end
        if ~isempty(nid)
            Header.excluded.delay = [Trials(nid).id];
        end
        Trials = Trials(id);
        durs = durs(id);
    end
    Expt.Trials = Trials;
    Expt.Header = Header;
    Expt.Header.trange(1) = AllExpts(nx).start;
    if isfield(AllExpts,'end')  && ~isempty(AllExpts(nx).end)
% If reaal End of expt is long time after last trial, probabaly means user did
% EndExpt. Still want to record real time of endexpt
        if length(AllExpts) < nexpts || Expt.Trials(end).End(end)+10000 > AllExpts(nx).end
            AllExpts(nx).end = Expt.Trials(end).End(end)+10000;
        else
            AllExpts(nx).end = AllExpts(nx).end+5000;
        end
    else %%online, unfinished, edpt
        if nx > 1  
            fprintf('Expt %d No end. Adding.\n');
        end
        AllExpts(nx).end = Expt.Trials(end).End(end)+10000;
    end
    if isfield(AllExpts,'end')  && ~isempty(AllExpts(nx).end)
        Expt.Header.trange(2) = AllExpts(nx).end;
    else
        Expt.Header.trange(2) = 0;
    end



    if isfield(Idx,'DigMark')
        t = Expt.Header.trange./10000;
        id = find(Idx.DigMark.times > t(1) & Idx.DigMark.times < t(2));
        if length(id)
            Expt.DigMark.times = Idx.DigMark.times(id);
            Expt.DigMark.codes = Idx.DigMark.codes(id,1);
        end
    end
    if psychtrial > nt/2 && sum(ismember(fields(Expt.Trials),{'sq' 'vs' 'sa' }))
        Expt.Header.psych = 1;
    else
        Expt.Header.psych = 0;
    end
    %if Expt.Stimvals.x2 is set, it means it was manually set via AddTxt
    if ~isfield(Expt.Stimvals,'x2') || Expt.Stimvals.x2 == 0
        xid = strcmp(checkoptions,'+x2');
        if optioncounts(xid) > nt/2
            Expt.Stimvals.x2 = 1;
        else
            Expt.Stimvals.x2 = 0;
        end
    end
    if seqtrial > nt/2
        Expt.Header.rc = 1;
    else
        Expt.Header.rc = 0;
    end
    Expt.Header.Options = '';
    for c = 1:length(checkoptions)
        if optioncounts(c) > nt/2
            Expt.Header.Options = [Expt.Header.Options checkoptions{c}];
        else
            Expt.Header.Options = [Expt.Header.Options strrep(checkoptions{c},'+','-')];
        end
    end
    if fastseq
        Expt = SetExptRC(Expt);
    elseif isfield(Expt.Trials,'rptframes')
        Expt.Trials = rmfields(Expt.Trials,{'mtFn' 'mtFl' 'mtFi'});
    end

    
    [a, expname, b, stimname] = Expt2Name(Expt,'addsuff');
    if isempty(expname)
        Expt.Header.expname = 'None';
    else
        Expt.Header.expname = a;
    end
    
    Expt.Header.Start = AllExpts(nx).start+timeoffset;
    Expt.Header.End = AllExpts(nx).end+timeoffset;
    Expt.Header.exptno = nexpts;    

    
    if min(durs) < 0
        id = find(durs < 0);
        nmiss = 0;
        xstr = '';
        if isfield(Idx,'frametimes') && isfield(Trials,'bstimes')
            gf = [0 0 ];
            for j = 1:length(id)
%Trials.Start is messed up by missing frametime, so use bstimes - the DIO time                
                nft = sum(Idx.frametimes > Trials(id(j)).bstimes(1) & Idx.frametimes < Trials(id(j)).End(end));
                if nft == 0 %no frametime signals in the trial. 
                    if gf(1) == 0
                        xid = find(Idx.frametimes < Trials(id(j)).bstimes(1));
                        if ~isempty(xid)
                            gf(1) = Idx.frametimes(xid(end));
                        end
                    end
                    xid = find(Idx.frametimes > Trials(id(j)).End(1));
                    gf(2) = Idx.frametimes(xid(1));
                    nmiss = nmiss+1;
                    cprintf('red','No Frame Signals for Trial Id %d - storge probably off\n',Trials(id(j)).id);
                end
            end
            if nmiss > length(id)/2;
                xstr = sprintf(': Storage off %.0f-%.0f',gf);
            end
        end
        err = sprintf('Expt %d at %.2f(id%d) has (%d) negative durations%s',nexpts,Expt.Trials(id(1)).Start(1),Expt.Trials(id(1)).id,length(id),xstr);
        Idx = AddError(Idx, err);
        id = find(durs > 0);
        Expt.Trials = Trials(id);
    end
    if ~isfield(Expt,'e1vals') & isfield(Expt.Trials,et)
        Expt.e1vals = unique(cat(1,Expt.Trials.(et)));
    end
 
    if ischar(Expt.Stimvals.ve)
        if strncmp(Expt.Stimvals.ve,'binoclean.',10)
            ve = sscanf(Expt.Stimvals.ve(11:end),'%f');detail
            Expt.Stimvals.ve = ve(1) .*10;
        end
    end
    if ~strcmp('e0',Expt.Stimvals.e2) && strcmp('dx',Expt.Stimvals.et) && Expt.Stimvals.ve < 10.1027
        Expt.Stimvals.x2 = 1;
    end
    
%replace off screen positiosn with 0, but only for 
%regular tuning curves
    if sum(strncmp(Expt.Stimvals.et,{'xo' 'yo' 'Op' 'Pp'},2)) && Expt.Header.rc == 0 && isfield(Expt.Trials,Expt.Stimvals.et)
        exv = [Expt.Trials.(Expt.Stimvals.et)];
        if size(exv,1) == 1
        id = find(abs(exv) > 35); %off screen
        for j = 1:length(id)
            Expt.Trials(id(j)).xo = Expt.Stimvals.rf(1);
            Expt.Trials(id(j)).yo = Expt.Stimvals.rf(2);
            Expt.Trials(id(j)).Op = 0;
            Expt.Trials(id(j)).Pp = 0;
            Expt.Trials(id(j)).st = 0;
        end
        end
    end
    
    if strcmp(Expt.Stimvals.et,'fx') && strcmp(Expt.Stimvals.e2,'fy')
        if isfield(Expt.Header,'loadname')
            emname = strrep(Expt.Header.loadname,'.mat','.eyecal.mat');
        else
            emname = strrep(Expt.Header.Name,'.mat','.eyecal.mat');
        end
        if ~exist(emname,'file')
            MakeEyeCalfile(Expt,emname);
        else
            a = load(emname);
            if isfield(a,'gains')
                fprintf('Eye gains RH %.2f LH %.2f RV %.2f LV %.2f\n',a.gains(1),a.gains(2),a.gains(3),a.gains(4));
            end
        end
    end
%     if strmatch(et,'Dc')
%         ors = cat(2,Expt.Trials.or);
%     end
    if isfield(Idx,'errexpt')
        id = find(Idx.errexpt == nexpts);
        if ~isempty(id)
            Expt.errs.errmsg = Idx.errs(id);
            if isfield(Idx,'errimtes')
                Expt.errs.t = Idx.errtimes(id);
            else
                Expt.errs.t = zeros(size(id));
            end
        end
    end
        Expts{nexpts} = Expt;
    nexpts = nexpts+1;
    end
%  max(spkids);
end
Expts = AddComments(Expts,Idx);
ExptState = rmfields(state,'toplevel');
Tidx.name = outname;
Tidx = CopyFields(Tidx, Idx, {'DataType' 'errs' 'state' 'ExptList'});
Tidx.state.savetime = now;
ExptState.nexpts = length(Expts);
save(outname,'Expts','ExptState','Tidx');
    
function  gains = MakeEyeCalfile(cExpt,emname);
 
if isfield(cExpt.Header,'loadname')
    emfile = strrep(cExpt.Header.loadname,'.mat','.em.mat');
else
    emfile = strrep(cExpt.Header.Name,'.mat','.em.mat');
end
if ~exist(emfile,'file')
    gains = NaN;
    return;
end
load(emfile);
if ~isfield(Expt.Trials,'fx') && isfield(Expt.Trials,'id')
    for j = 1:length(cExpt.Trials)
        id = find([Expt.Trials.id] == cExpt.Trials(j).id);
        Expt.Trials(id).fx = cExpt.Trials(j).fx;
        Expt.Trials(id).fy = cExpt.Trials(j).fy;
    end
end

if isfield(Expt.Trials,'fx') && isfield(Expt.Header,'emtimes')
    [gains, positions] = EyeCal(Expt,'ids',[cExpt.Trials.id]);
    save(emname,'gains','positions');
end

function d = CellCreationDate(Text)

    did = find(strncmp('uf',Text.text,2));
if isempty(did) %online file
    did = find(strncmp('bt',Text.text,2));
end
d = 0;
if length(did)
    for j = 1:length(did)
        ds = Text.text{did(j)};
        dsid = strfind(ds,'Creat');
        if length(dsid)
            d = datenum(ds(dsid(1)+8:end));
            break;
        end
    end
end
if d == 0 %still not found
    did = find(strncmp('vebinoc',Text.text,7));
    if ~isempty(did)
        ds = Text.text{did(j)};
        did = strfind(ds,' ');
         if length(did)
            d = datenum(ds(did(1)+1:end));
        end
    end
end

function d = CreationDate(Text)

if iscell(Text.text)
    d = CellCreationDate(Text);
    return;
end
did = strmatch('uf',Text.text);
if isempty(did) %online file
    did = strmatch('bt',Text.text);
end
d = 0;
if length(did)
    for j = 1:length(did)
        ds = Text.text(did(j),:);
        dsid = strfind(ds,'Creat');
        if length(dsid)
            d = datenum(ds(dsid(1)+8:end));
            break;
        end
    end
end
if d == 0 %still not found
    did = strmatch('vebinoc',Text.text);
    if ~isempty(did)
        ds = Text.text(did(1),:);
        did = strfind(ds,' ');
         if length(did)
            d = datenum(ds(did(1)+1:end));
        end
    end
end


function [rfstr, rf] = MkUfl(name, Text, varargin)
%first make .ufl file with rf boxes, so that can build pen maps the
%old way
overwrite = 0;
j = 1;
rfstr = [];
rf = [];
while j <= length(varargin)
    if strncmpi(varargin{j},'overwrite',5)
        overwrite = 1;
    end
    j = j+1;
end
ufl = strrep(name,'.mat','.ufl');
if iscell(Text.text)
    rid = find(strncmp('cm=rf',Text.text,5));
else
    rid = strmatch('cm=rf',Text.text);
end



AddTxtFile = strrep(name,'.mat','Add.txt');
if ~exist(AddTxtFile)
    AddTxtFile= regexprep(name,'\.[0-9]*.mat','Add.txt');
end
fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%d %s','delimiter','\n');
    fclose(fid);
    id = find(a{1} == -1);
    rfstrs = a{2}(id);
else
    rfstrs = {};
end

if isempty(rid) && isempty(rfstrs)
    return;
end
%a = textscan(Text.text(id,:),'cm=rf%f,%f:%fx%f,%fdeg pe%d %f %f%*s');
% trailing spaces seem to mess this up. text(id,1:65) works for most line
% but still barfs if a line is the wrong length
if isempty(rfstrs)
for j = 1:length(rid)
    if iscell(Text.text)
    a = sscanf(Text.text{rid(j)}','cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
    else
    a = sscanf(Text.text(rid(j),:)','cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
    end
    rfs(j,1:length(a)) = a;
end
else
    for j = 1:length(rfstrs)
        a = sscanf(rfstrs{j},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
        rfs(j,1:length(a)) = a;
    end
end
% find lines suggesting RF was changed after a quantitative measure
if iscell(Text.text)
    oid = find(strncmp('RO',Text.text,2))';
    pid = find(strncmp('RP',Text.text,2))';
    sid = [oid pid find(strncmp('RO',Text.text,2))'];
else
    oid = strmatch('RO',Text.text);
    pid = strmatch('RP',Text.text);
    sid = [oid(:)' pid(:)'];
end
if length(sid)
    id = find(rid > max(sid))
end
for j = 1:size(rfs,2)
    rf(j) = mode(rfs(:,j));
end
if size(rfs,2) < 10
  rfstr = 'Missing RF data';
else
    rfstr = sprintf('cm=rf%.2f,%.2f:%.2fx%.2f,%.0fdeg pe%.0f %.1f,%.1f fx=%.2f,fy=%.2f\n',...
        rf(1),rf(2),rf(3),rf(4),rf(5),mode(rfs(:,6)),...
        mode(rfs(:,7)),mode(rfs(:,8)),mode(rfs(:,9)),mode(rfs(:,10)));
end




if exist(ufl,'file') & ~overwrite
    return;
end


d = CreationDate(Text);
ds = [];
if d > 0
    ds = ['Created: ' datestr(d,'mm/dd/yyyy')];
end
of = fopen(ufl,'w');
if of > 0 
    fprintf(of,'%s\n',rfstr);
    for j = 1:length(rid)
%        fprintf(of,'%s\n',Text.text(id(j),:));
    end
    if ~isempty(ds)
        fprintf(of,'%s\n',ds);
    end
    fclose(of);
else
    questdlg(sprintf('Can''t Write %s',ufl),'test','OK','OK');
end

function Idx = ShowErrors(Idx, Spkerrs)

if isempty(Spkerrs)
    return;
end
errstr = {Spkerrs.s};
strs = unique(errstr);
for j = 1:length(strs)
    id = find(strcmp(strs{j},errstr));
    pstr = sprintf('%d,',unique([Spkerrs(id).p]));
    pstr(end) = ':';
    strs{j}(strs{j} < 12) = []; %remove non-printing
    if strfind(strs{j},':Not Applying')
        show = '-silent';
    else
        show = '-noshow';
    end
    Idx = AddError(Idx,show,'P%s %s',pstr,deblank(strs{j}));
end

function Idx = AddError(err, varargin)
%Idx = AddError(Idx, varargin)
% or Idx = AddError(err, Idx, show) old style
silent = 0;
if ischar(err) %% old style
    Idx = varargin{1};
    show = varargin{2};
else  
    Idx = err;
    if isfield(Idx,'state')
        show = Idx.state.showerrs;
    else
        show = 0;
    end
    if sum(strncmp(varargin{1},{'-noshow' '-silent'},6)) %log error, but don't pop warning
        show = 0;
        if strcmp(varargin{1},'-silent')
            silent = 1; %don't print these
        end
        varargin = varargin(2:end);
    elseif strncmp(varargin{1},'-nonmodal',6) %don't pop modal warning
        show = 2;
        varargin = varargin(2:end);
    end
    err = sprintf(varargin{:});    
end
        
if ~isfield(Idx,'errs') | sum(strncmp(err,Idx.errs,length(err))) == 0 
    if show ==2
        h = msgbox(err,'APlaySpkFile Error!!');
    elseif show == 1
        h = msgbox(err,'APlaySpkFile Error!!','modal');
    else
        h = 0;
    end
    if myhandle(h)
        setappdata(h,'ParentFigure','APlaySpkFile');
    end
    if silent == 0
        mycprintf('errors',[strrep(err,'\','/') '\n']);
    end
    if ~isfield(Idx,'errs')
        Idx.errs{1} = err;
    else
        Idx.errs = {Idx.errs{:} err};
    end
    nerr = length(Idx.errs);
    Idx.errdata(nerr).s = err;
    if isfield(Idx,'t')
        Idx.errtimes(nerr) = Idx.t;
        Idx.errdata(nerr).t = Idx.t;
    end
    if isfield(Idx,'exptno')
        Idx.errexpt(nerr) = Idx.exptno;
        Idx.errdata(nerr).exptno = Idx.exptno;
    end
    Idx.errdata(nerr).time = now;
    Idx.errdata(nerr).progname = 'APlaySpkFile';
    if isfield(Idx,'newerrs')
        Idx.newerrs = Idx.newerrs+1;
    end
else
    mycprintf('errors',[err '\n']);
end

function s = FindMissingTimes(Events, Text, bstimes, estimes,bsid,fsid,ts,te)

s = '';
stimid = 0;
bid = strmatch('bss',Text.text);
bsstimes = Text.times(bid);
tid = find(bsstimes >  ts & bsstimes < te);
eid = find(bstimes > ts & bstimes < te);
nl = 1;
if length(bstimes) >= length(bsstimes)
x = bsstimes-bstimes(1:length(bsstimes));
bdiff = diff(bstimes(bsid));
ediff = diff(Events.times(fsid));
subplot(2,1,2);
hold off; 
plot(bdiff);
hold on;
plot(ediff);
[n, short] = min([length(bdiff) length(ediff)]);
jdiff = bdiff(1:n)-ediff(1:n);
emax = min([max(abs(jdiff))/2 0.2e5]);
id = find(abs(jdiff) >emax,1);
if jdiff(id) > 0 %bstime has big jump - its missing one
    s = sprintf('Missing DIO On at about %.4f id%d',Events.times(fsid(id))./10000,stimid);
    fprintf('Looks like a missing bstime at %.4f\n',Events.times(fsid(id))./10000);
else
    xs = (bstimes(bsid(1:length(fsid)))-Events.times(fsid))./10000;
    tx(1:2) = bstimes(bsid(id:id+1))./10000;
    tx(3:4) = Events.times(fsid(id:id+1))./10000;
    fprintf('Looks like a missing Frame signal %.4f. Next is %.4f\n',tx(2),tx(4));
    idid = find(strncmp('id',Text.text,2) & Text.times./1000 > tx(2)-0.1,1);
    stimid = sscanf(Text.text{idid},'id%d');
    fprintf('Previous pair %.4f and %.4f = %.4f diff\n',tx(1),tx(3),diff(tx([1 3])));
    if id >3 && id < length(fsid) - 2
        idmax = min([length(bsid) length(fsid) id+4]);
        fprintf('Frame :')
        for j = id-2:idmax
            fprintf('%.3f  ',Events.times(fsid(j))./10000);
        end
        fprintf('\n');
        fprintf('StimOn:')
        for j = id-2:idmax
            fprintf('%.3f  ',bstimes(bsid(j))./10000);
        end
        fprintf('\n');
        fprintf('StimOff:')
        for j = id-2:idmax
            fprintf('%.3f  ',estimes(bsid(j))./10000);
        end
        fprintf('\n');
        
    end
    a = corrcoef(bdiff(1:id),ediff(1:id));
    b = corrcoef(bdiff(id+2:length(ediff)+1),ediff(id+1:end));
    if length(bdiff) > length(ediff)+1
        c = corrcoef(bdiff(id+3:length(ediff)+2),ediff(id+1:end));
        if c(1,2) > b(1,2) && c(1,2) > 0.9
           s = sprintf('Need to add two codes 5  at %.4f StimON dur was %.4f id%d\n',tx(2),(estimes(id)./10000)-tx(2),stimid);
        end
    end
    if abs(diff(tx([2 4]))) > 0.5 && abs(diff(tx([1 3]))) < 0.1
        id = find(estimes./10000 > tx(2),1);
        fprintf('interval xcorr before/after with fix: %.3f,%.3f\n',a(1,2),b(1,2));
        s = sprintf('Should be safe to remove StimOn (see Spike2.FileHelp) at %.4f StimON dur was %.4f id%d\n',tx(2),(estimes(id)./10000)-tx(2),stimid);
    end
end

subplot(2,1,1);
hold off;
plot(bsstimes,x,'o');
hold on;
labels{nl} = 'txt - ev (all)';
nl = nl+1;
end
if length(tid) <= length(eid) 
%seem to be missing serial line Frame signal evenr and or bss text    
    x = bsstimes(tid)-bstimes(eid(1:length(tid)));
    plot(bsstimes(tid),x,'ro');
    labels{nl} = 'txt - ev';
    nl = nl+1;
    xid = find(diff(x) > 10000);
    if ~isempty(xid)
        t = bstimes(eid(xid(1)+1));
        fprintf('Unmatched StimOn at %.3f\n',t./10000);
    end
else
    x = bsstimes(tid(1:length(eid)))-bstimes(eid);
    plot(bsstimes(tid(1:length(eid))),x,'ro');
    labels{nl} = 'txt - ev (in expt)';
    nl = nl+1;
end
plot([ts ts],get(gca,'ylim'));
plot([te te],get(gca,'ylim'));
if length(fsid)  > length(bsid)
for j = 1:length(bsid)
    adiffs(j) = bstimes(bsid(j)) - Events.times(fsid(j));
    bdiffs(j) = bstimes(bsid(j)) - Events.times(fsid(j+1));
end
cdiff = abs(adiffs) > abs(bdiffs);
%if adiffs is small at first, then it changes so that bdiffs is small
%thats the missing trial
id = 1+find(diff(cdiff) > 0);
if length(id) == 1
    cprintf('red','MissingStimon For Trial %d at %.3f\n',id, Events.times(fsid(id)));
    cprintf('red','Create .Stimon File to fix\n');
end
plot(adiffs,bdiffs);
end
if length(bsid) > length(fsid)
    for j = 1:length(fsid)
        adiffs(j) = bstimes(bsid(j)) - Events.times(fsid(j));
        bdiffs(j) = bstimes(bsid(j+1)) - Events.times(fsid(j));
    end
    cdiff = abs(adiffs) > abs(bdiffs);
end
legend(labels);

function name = BuildName(name)
if isempty([strfind(name,'/') strfind(name,'\')])
    name = [pwd '/' name];
end
id = strfind(name,':');
if id
    name = name(id(1)+1:end);
else
    name = name;
end
if isunix
    name = strrep(name,'\','/');
end


