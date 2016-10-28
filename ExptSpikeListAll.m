function espk = ExptSpikeListAll(DATA, eid, stimes)
    t(1) = DATA.Expts{eid}.Trials(1).Start(1)-DATA.state.preperiod;
    t(2) = DATA.Expts{eid}.Trials(end).End(end)+DATA.state.postperiod;
    if DATA.state.online && (max(stimes) - t(2) < 100000 || sum(stimes>t(2) < 1000))%<10 sec after last trial = use all
        espk = find(stimes > t(1));
    else
        espk = find(stimes > t(1) & stimes < t(2));
    end