function C = CheckClusterFields(C, varargin)nprobes = 24; %defaultDATA = [];j = 1; while j <= length(varargin)    if isfield(varargin{j},'autocutmode')        DATA = varargin{j};    end    j = j+1;endif iscell(C)    for j = 1:length(C)        C{j} = AllV.CheckClusterFields(C{j}, varargin{:});    end    return;endif isfield(C,'MeanSpike') && size(C.MeanSpike.ms,1) > 4    nprobes = size(C.MeanSpike.ms,1);end    if ~isfield(C,'auto')        C.auto = 0;    end    if ~isfield(C,'clusterprog')        C.clusterprog  = 'AllVPcs 1.00';    end    if ~isfield(C,'manual')        C.manual = 0;    end    if ~isfield(C,'probe')        if isfield(C,'chspk')            C.probe = c.chspk(1);        else            C.probe = 0;        end    end    if ~isfield(C,'chspk')        if isfield(DATA,'chspk') && DATA.probe == C.probe            C.chspk = DATA.chspk;        else            C.chspk = C.probe + [-1:1];            C.chspk = C.chspk(C.chspk > 0 & C.chspk <= nprobes);        end    end    if ~isfield(C,'triggerchan')        C.triggerchan = C.probe(1);    end    if ~isfield(C,'tsmooth')        C.tsmooth = 0;;    end    