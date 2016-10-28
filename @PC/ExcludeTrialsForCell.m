function DATA = ExcludeTrialsForCell(DATA, probe, cluster, varargin)        useall = 0;    j = 1;    while j <= length(varargin)        if strncmpi(varargin{j},'reset',5)            useall = 1;        end       j = j+1;    end    Clusters = getappdata(DATA.toplevel,'Clusters');    e = DATA.currentpoint(1);    C = Clusters{e}{probe};    if probe == 0        it = findobj('Tag','CellNumberId');        cellid = get(it,'value');    else        cellid = DATA.CellList(e, probe, cluster);    end%    if cellid == 0%        return;%    end           F = FindFig(DATA.tag.spikes);    if length(F) == 1        it = findobj(F,'Tag','ChooseTrial');    else        it = findobj(DATA.spoolfig,'Tag','TrialList');    end    trials = get(it,'value');    oldt = PC.FindExcludedTrials(DATA, e, probe, cluster,C);    trials = unique([trials oldt]);    if useall        DATA.CellDetails.excludetrials{ e, probe, cluster} = [];    else        DATA.CellDetails.excludetrials{e, probe, cluster} = DATA.trialids{e}(trials);    end    fprintf('Cell%d (E%dP%dCluster%d) Excluding: %s\n',cellid,e,probe,cluster,sprintf(' %d',trials));set(DATA.toplevel,'UserData',DATA);DATA = PC.SaveCellList(DATA);if length(F) == 1    PC.SetTrialList(DATA,C,DATA.currenttrial);else    PC.SpoolSpikes(DATA.spoolfig,'excludelist', trials);end