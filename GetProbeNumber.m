function p = GetProbeNumber(DATA, varargin)
%p = GetProbeNumber(DATA) extracts probe # from struct/string

if ischar(DATA);
    p = GetProbeFromName(DATA, varargin{:});
elseif iscell(DATA)
    for j = 1:length(DATA)
        p(j) = GetProbeNumber(DATA{j});
    end
elseif isfield(DATA,'trueprobe') && DATA.trueprobe > 0
    p = DATA.trueprobe;
elseif isfield(DATA,'probelist')
    if DATA.probe(1) > length(DATA.probelist)
        p = DATA.probe(1);
    elseif DATA.probe(1) == 0
        p = DATA.probelist(1);
    else
        p = DATA.probelist(DATA.probe(1));
    end
elseif isfield(DATA,'probe')
    p = DATA.probe(1);
elseif isfield(DATA,'Header')
    if isfield(DATA.Header,'probe')
        p = DATA.Header.probe;
    else
        p = NaN;
    end
else
    p = NaN;
end
    