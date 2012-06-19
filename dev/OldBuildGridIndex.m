function idx = BuildGridIndex(name, Expts, varargin)
%idx = BuildGridIndex(name, Expts, ...)
%Build an index of which neV files match Expt .mat files
%name is a filename or direcotry where the .mat file lives (used to
%construct path for where .nev files live
% Expts is a cell array of Expts, as returned by APlaySpkFile
    idx = [];
   datdir = 'F:/Utah/dufus/';
   reindex = 0;
   plotexpts = [];
   checktag = 'TrialCheck';
   preperiod = 5000;
   postperiod = 5000;
   maxgap = 10000; %default is to chop if gap > 1 sec
   DigMark = [];
   
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'digmark',5)
        j = j+1;
        DigMark = varargin{j};
    elseif strncmpi(varargin{j},'plotexpts',5)
        j = j+1;
        plotexpts = varargin{j};
    elseif strncmpi(varargin{j},'reindex',5)
        reindex = 1;
    end
    j = j+1;
end
   
   if isdir(name)
       datdir = name;
   else
       datdir = fileparts(name);
   end
   plotsummary = 2;

   if isempty(strfind(path,'BlackRock'))
       path(path,'/bgc/matlab/BlackRock');
   end
   idxfile = [datdir '/FileIdx.mat'];
   
   if exist(idxfile,'file') & ~reindex
       a = load(idxfile);
       idx = a.idx;
       idx.datdir = datdir;
       nidx = length(idx.names);
       return;
   else
       nidx = 0;
   end
   
   
   if isempty(Expts)
       idx = BuildMatFiles(datdir);
       return;
   end
   
   
   offtimes = [];
   ontimes = [];
   sampleoffs = [];
   sampleons = [];
   offids = [];
   ncut = 0;
% First fix bug with early online files where DigMark got inverted somehow
%should be able to remove these lines one day.....
   for j = 1:length(Expts)
       if length(Expts{j}.gridstoreon) ==1 && length(Expts{j}.DigMark.times) == 1 ...
               && Expts{j}.gridstoreon > Expts{j}.gridstoreoff
           a = Expts{j}.gridstoreoff;
           Expts{j}.gridstoreoff = Expts{j}.gridstoreon;
           Expts{j}.gridstoreon = a;
           if Expts{j}.DigMark.codes == 1
               Expts{j}.DigMark.codes  = 2;
           end
       end
   end
   
   for j = 1:length(Expts)
       starts(j) = Expts{j}.Header.CreationDate + Expts{j}.Header.Start./(10000 * 60 *60 *24);
       ebstimes{j} = Expts{j}.bstimes;
       id = find(Expts{j}.DigMark.codes ==1);
       if isempty(id)
       ontimes = [ontimes Expts{j}.gridstoreon./10000];
       else
           %are digmark times in sec or 1/10 msec? Online may differ from
           %final
       ontimes = [ontimes Expts{j}.DigMark.times(id(1))];
       end
       if id == 1
           sampleon = Expts{j}.Header.CreationDate + Expts{j}.DigMark.times(id)'./(60 *60 *24);
       elseif isempty(id)
           sampleon = Expts{j}.Header.CreationDate + Expts{j}.gridstoreon./(1000 * 60 *60 *24);
       else
           sampleon = Expts{j}.Header.CreationDate + Expts{j}.DigMark.times(id)'./(60 *60 *24);
       end
       sampleons = [sampleons sampleon];

       id = find(Expts{j}.DigMark.codes ==2);
       if isempty(id)
           fprintf('Expt %d Missing Sample off marker\n',j);
           sampleoff = Expts{j}.Header.CreationDate + Expts{j}.gridstoreoff(1)./(10000 * 60 *60 *24);
           offtimes = [offtimes Expts{j}.gridstoreoff(1)./10000]
       else
           sampleoff = Expts{j}.Header.CreationDate + Expts{j}.DigMark.times(id)'./(60 *60 *24);
           offtimes = [offtimes Expts{j}.DigMark.times(id)'./10000];
       end
       sampleoffs = [sampleoffs sampleoff];
       offids = [offids ones(size(sampleoff)).*j];
%       eestimes{j} = Expts{j}.estimes;
       for t = 2:length(Expts{j}.Trials)
           gaps(t) = Expts{j}.Trials(t).Start(1)-Expts{j}.Trials(t-1).End(end);
           if  gaps(t) > maxgap
            ncut = ncut+1;   
            cuts(ncut,1) = Expts{j}.Trials(t-1).End(end)+postperiod;
            cuts(ncut,2) = Expts{j}.Trials(t).Start(1)-preperiod;
           end
       end
   end
   ends(1:j-1) = starts(2:end);
   ends(j) = starts(j) + (max(ebstimes{j})./(10000 * 60 *60 *24));
   
   
%First Build a list of Nev files and their times
%This works for multiple files per expt. 
   d = dir([datdir]);
   filenames = {d.name};
   newf = 0;
   nnev = 0;
   for j = 1:length(d)
       if strfind(d(j).name,'.nev') %this has starttime and Dig Events
           nevfile = [datdir '/' d(j).name];
           matfile = strrep(d(j).name,'.nev','.mat');
           mid = strmatch(matfile,filenames);
           if length(mid)
               agediff = d(j).datenum > d(mid).datenum;
           else
               agediff = 1;
           end
           if reindex || isempty(idx) || isempty(strmatch(d(j).name,idx.names)) || ...
                   agediff > 0
               if  agediff > 0
                   nev = openNEV('read','nomat','noparse','nowarning',nevfile);
               else
                   nev = openNEV('read','noparse','nowarning','nowaves',nevfile);
               end
               eb = int16(bitand(3,nev.Data.SerialDigitalIO.UnparsedData));
               onoff = diff(int16(bitand(4,nev.Data.SerialDigitalIO.UnparsedData)));
               nnev =nnev+1;
               id = find(onoff > 0);
               bstimes = nev.Data.SerialDigitalIO.TimeStampSec(id+1).*10000;
               id = find(onoff < 0);
               estimes = nev.Data.SerialDigitalIO.TimeStampSec(id+1).*10000;
               ns = length(bstimes);

               ts = nev.MetaTags.DateTimeRaw;
               tstart  = datenum(ts(1),ts(2),ts(4),ts(5),ts(6),ts(7));
%When storage is turned off the lowest bit is dropped 1ms boefore bit2, 
% so lowest two bits == 2. Only other time this happens is at the stare
% when bit 1 is toggled to mark times.
               if length(nev.Data.SerialDigitalIO.UnparsedData) == 0 || ...
                       bitand(3, nev.Data.SerialDigitalIO.UnparsedData(end)) > 0
                   fprintf('%s Missing Final Off event',nevfile);
                   [a,c] = min(abs(tstart - sampleons));
                   b = offids(c);
                   cpuclockdiff = (tstart-sampleons(c)).*(60*60*24);
                   tstop = tstart;
                   nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(end);
                   tdiff = ontimes(c)-nevoff; %difference in sec
                   cpuclockdiff = (tstop-sampleoffs(c)).*(60*60*24);

               else
                   id = find(bitand(3,nev.Data.SerialDigitalIO.UnparsedData) ==2);
                   nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(id(end));
                   if nevoff < 1
                       fprintf('Nominally short file %.2f\n', nevoff);
                       nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(end);
                       if isfield(nev.Data.Spikes,'TimeStamp')
                           lastspk = double(nev.Data.Spikes.TimeStamp(end))./nev.MetaTags.TimeRes;
                           nevoff = max([nevoff lastspk]);
                       end
                       t = tstart+nevoff./(60 * 60 * 24);
                       xid = find(sampleoffs) > t;
                       fprintf('Missing ~ %.2f sec of data\n',sampleoffs(id(1))-1);
                   end
                   t = tstart+nevoff./(60 * 60 * 24);
                   tstop = t;
                   [a,c] = min(abs(t - sampleoffs));
                   b = offids(c);
                   tdiff = offtimes(c)-nevoff; %difference in sec
                   cpuclockdiff = (tstop-sampleoffs(c)).*(60*60*24);
               end
               requesttime = sampleons(c) * 10000;
               nbs = length(ebstimes{b});
               E = Expts{b};
               nevstarts(nnev) = tstart;
               nevfiles{nnev} = d(j).name;
               if ns && nbs
                   diffs = ebstimes{b} - (bstimes(1)+tdiff*10000);
                   [a,bsoff] = min(abs(diffs));
               end

               if ns && length(estimes) && nbs && ...
                       abs(diffs(bsoff)) < 10000 && abs(cpuclockdiff) < 4
                   nidx = nidx+1;
                   newf = newf+1;
                   idx.names{nidx} = d(j).name;

                   idx.expt(nidx) = b;
                   idx.tdiff(nidx) = (tstop-sampleoffs(c)).*(60*60*24);

                   
                   extimes = ebstimes{b}(bsoff:bsoff+length(bstimes)-1);
                   xc = corrcoef(diff(bstimes),diff(extimes));
%toff is the timestamp in spike2 associated with t = 0 in the Cerebrus file
                   idx.toff(nidx) = ebstimes{b}(bsoff)-bstimes(1);
                   [a, subid] = min(abs(idx.toff(nidx)-requesttime));
%delay is best estimate of the delay between the request from spike2 and the file opending
%in Cerebrus. Based where t=0 falls in spike2, relative to requesttime (dig
%marker in spike2). ddelay (below) is time of first recorded event. Can be
%later if all timing pulses are missed
                   idx.delay(nidx) = idx.toff(nidx)-requesttime(subid);
                   idx.firstbs(nidx) = bsoff;
                   idx.start(nidx) = tstart + bstimes(1)./(10000 * 60 * 24); %datenum
                   idx.end(nidx) = tstart + estimes(end)./ (10000 * 60 * 24);
                   idx.nt(nidx) = length(estimes);
                   idx.bstimes{nidx} = bstimes;
                   idx.evtimes{nidx} = nev.Data.SerialDigitalIO.TimeStampSec .* 10000;
                   idx.digin{nidx} = bitand(7,nev.Data.SerialDigitalIO.UnparsedData);
                   idx.digdt(nidx) = diffs(bsoff);
                   cid = find(cuts(:,1) > idx.toff(nidx) & cuts(:,1) < idx.toff(nidx)+nevoff*10000);
                   idx.cuts{nidx} = cuts(cid,:);
                   
               else
                   nidx = nidx+1;
                   idx.expt(nidx) = 0;
                   idx.tdiff(nidx) = NaN;
                   if ns
                   idx.digdt(nidx) = diffs(bsoff);
                   else
                   idx.digdt(nidx) = NaN;
                   end
               end
               if length(nev.Data.SerialDigitalIO.TimeStampSec)
                   idx.ddelay(nidx) = nev.Data.SerialDigitalIO.TimeStampSec(1).*10000;
               else
                   idx.ddelay(nidx) = NaN;
               end
               if length(nev.Data.SerialDigitalIO.TimeStampSec)
                   idx.ddelay(nidx) = nev.Data.SerialDigitalIO.TimeStampSec(1).*10000;
               else
                   idx.ddelay(nidx) = NaN;
               end
               if c > 1
                   idx.offinterval(nidx) = sampleons(c)-offtimes(c-1);
               end
               idx.starttime(nidx) = tstart;
               idx.stoptime(nidx) = tstop;
               idx.names{nidx} = d(j).name;
               if length(nev.Data.SerialDigitalIO.UnparsedData) 
                   idx.lastdio(nidx) = bitand(7,nev.Data.SerialDigitalIO.UnparsedData(end));
               end
           end
       end
   end
   if isempty(idx)
       fprintf('Missing NeV Data Filesin %s\n',name);
       return;
   end

   if newf
       save(idxfile,'idx');
   end
idx.datdir = datdir;

if length(plotexpts)
    GetFigure(checktag);
    hold off;
    for j = 1:length(plotexpts)
        PlotTrialMarks(idx, Expts, plotexpts(j));
        hold on;
    end
end

   
   
function PlotTrialMarks(idx, Expts,  exptno)
   
        plotsummary = 2;
   id = find(idx.expt == exptno);
   bid = id;
   for j = 1:length(id)
       nfiles{j} = [idx.datdir '/' idx.names{id(j)}];
   end
   if plotsummary == 2 %plot trial start/end
       b = exptno;
       T = Expts{b}.Trials;
       ebstimes = Expts{b}.bstimes;
    for j = 1:length(ebstimes)
%        plot([ebstimes(j) ebstimes(j) eestimes{b}(j) eestimes{b}(j)],[0 1 1 0],'r-');
        plot([ebstimes(j) ebstimes(j)],[0 1],'r-');
        hold on;
    end
    for j = 1:length(T)
        plot([T(j).Start(1) T(j).Start(1)],[-1 -2],'g-');
    end
    for j = 1:length(id)
        bs = idx.bstimes{id(j)} + idx.toff(id(j));
        for k = 1:length(bs)
            plot([bs(k) bs(k)],[0 -1 ],'-');
        end
        text(mean(bs), -1,idx.names{id(j)}(8:10));
        for k = 1:length(idx.evtimes{id(j)})
            t = idx.evtimes{id(j)}(k) + idx.toff(id(j));
            if bitand(1,idx.digin{id(j)}(k))
                plot([t t],[-0.5 -1.5],'g');
            else
                plot([t t],[-0.5 -1.5],'r');
            end
        end
    end
%    aid = find(idx.starttime > Expts{b}.Headerstarts(b) & idx.starttime < ends(b));
    aid = find(idx.toff > Expts{b}.Header.Start & idx.toff < Expts{b}.Header.End);
    for j = 1:length(aid)
        t = (idx.starttime(aid(j))-Expts{b}.Header.CreationDate) * (60 * 60 * 24 * 10000);
        plot([t t],[-2 -3],'r');
        dy = mod(j,5)/5;
        text(t, -2 - dy,idx.names{aid(j)}(8:10));
    end

    for j = 1:length(aid)
        t = idx.toff(aid(j));
        plot([t t],[-2 -3],'b');
        plot(t+idx.ddelay(aid(j)),-2.5,'x');
        dy = mod(j,5)/5;
        text(t, -2 - dy,idx.names{aid(j)}(8:10));
    end
    E = Expts{b};
    requesttime = [];
    if isfield(E,'DigMark')
        for j = 1:length(E.DigMark.times)
            t = E.DigMark.times(j).*10000;
            if bitand(1,E.DigMark.codes(j))
                plot([t t],[-4 -5],'r');
                requesttime = [requesttime t];
            else
                plot([t t],[-4 -5],'g');
            end
        end
    end
    t = E.gridstoreon(1);
    plot([t t],[-3 -4],'r');
    for j = 1:length(E.gridstoreoff)
        id = find(E.gridstoreon < E.gridstoreoff(j));
        if length(id)
        t = E.gridstoreon(id(end));
        plot([t t],[-3 -4],'r');
        end
    end
   end
   
   
   

                

function idx = BuildMatFiles(datdir)

idx = [];
reindex = 0;
nf = 0;

d = dir(datdir);
filenames = {d.name};
newf = 0;
for j = 1:length(d)
    if strfind(d(j).name,'.nev') %this has starttime and Dig Events
        nevfile = [datdir '/' d(j).name];
        matfile = strrep(d(j).name,'.nev','.mat');
        mid = strmatch(matfile,filenames);
        if length(mid)
            agediff = d(j).datenum > d(mid).datenum;
        else
            agediff = 1;
        end
        if reindex || isempty(idx) || isempty(strmatch(d(j).name,idx.names)) || ...
                agediff > 0
            nf = nf+1;
            if  agediff > 0
                nev = openNEV('read','nomat','noparse','nowarning',nevfile);
            else
                nev = openNEV('read','noparse','nowarning','nowaves', nevfile);
            end
            idx.names{nf} = d(j).name;
            idx.expt(nf) = nf;
            if ~isempty(nev.Data.Spikes.Electrode)
            idx.nprobes(nf) = max(nev.Data.Spikes.Electrode);
            end
            idx.toff(nf) = 0;
        end
    end
end
idx.datdir = datdir;
