function E = CountSpikes(Expt, C, clnum,xcl, varargin)%E = PC.CountSpikes(Expt, C, clnum,xcl)%finds spikes in C that match trials in Expt, and%adds them to the Expt.Trials strucutre.%clnum 1 = cluster 1, i.e. C.clst == 2strargs = cell2cellstr(varargin);fillmu = sum(strcmp('fillmu',strargs));preperiod = 1000;postperiod = 1000;        latency = 500;                id = find(C.clst == clnum+1);        if isfield(C,'t')            t = (C.t(id) .*10000);        else            t = C.times .* 1000;        end        if fillmu            muid = find(C.clst ~= clnum+1);            mutimes = (C.t(muid) .*10000);        end        for j = 1:length(Expt.Trials)            starts(j) =  Expt.Trials(j).Start(1);            ends(j) =  Expt.Trials(j).End(end);        end        useall = 1; % make user exclude thses trials. Otherwise not recorded        if isempty(id) || useall            npre = 0;            npost = 0;        else        npre = find(starts < t(1));        if npre > 10            id = find(starts < t(1)+10000);            Expt.Trials = Expt.Trials(id);        end        npost = find(starts > t(end));        if npost > 1            id = find(starts < t(end));            Expt.Trials = Expt.Trials(id);        end        end                        for j = 1:length(Expt.Trials)            T = Expt.Trials(j);            id = find(t > T.Start(1)-preperiod & t <= T.End(end)+postperiod);            Expt.Trials(j).Spikes = round([t(id)-T.Start(1)]');            Expt.Trials(j).count = sum(t(id) > T.Start(1)+latency & t(id) < T.End(end)+latency);            if size(t,2) == 1                Expt.Trials(j).Spikes = Expt.Trials(j).Spikes';            end            if fillmu                id = find(mutimes > T.Start(1)-preperiod & mutimes <= T.End(end)+postperiod);                Expt.Trials(j).OSpikes = round([mutimes(id)-T.Start(1)]');                Expt.Trials(j).Ocodes = C.clst(muid(id))-1;            end        end        E = Expt;        [a,b] = setdiff([Expt.Trials.id],xcl);        E.Trials = Expt.Trials(b);        