function [DATA, Vall, ispk, newdata] = SetupVall(DATA,Vall, ispk, newdata)if isfield(Vall, 'Spikes')    return;endif isfield(Vall,'t')    vt = Vall.t;elseif isfield(DATA,'blklen')    first = 1;    vt(sum(length(DATA.blklen))) = 0;    for j = 1:length(DATA.blklen)        last = first+DATA.blklen(j)-1;        vt(first:last) = DATA.blkstart(j)+[1:DATA.blklen(j)].*DATA.interval;        first = last+1;    end    if length(vt) > size(Vall.V,2)        PrintMsg(DATA.logname,'sum blklen > length(V)\n');        vt = vt(1:size(Vall.V,2));    end    if size(Vall.V,2) > length(vt)        vt(size(Vall.V,2)) = 0;    end    Vall.t = vt;else    vt= [1:size(Vall.V,2)] .* DATA.interval;    Vall.t = vt;endforceispk = 1; if isfield(DATA,'DataType')    if strncmp(DATA.DataType,'GridData',8)        DATA.probelist = 1:96;    endendif size(Vall.V,1) == 1 && ispk(1) > 0    newdata = 1;    if ~isempty(DATA.allints);        DATA.probelist = ispk; % for now        Vall.chspk = ispk;        DATA.trueprobe = ispk;        chspk = ispk;        ispk = 1;        DATA.chspk = chspk;        DATA.allnprobes = 64;        DATA.DataType = 'MicroWire';    elseif ispk > 1        DATA.probelist = ispk;        if ispk == Vall.chspk            ispk = 1;        elseif forceispk            Vall.chspk = ispk;            ispk = 1;        end    else        DATA.probelist = Vall.chspk;    endelse    DATA.probelist = 1:size(Vall.V,1);endif DATA.subtractmeanV    scale = (meanV*double(Vall.V'))./(meanV*meanV');    Vall.V = Vall.V - int16(meanV.*scale);endpres = {};if size(Vall.V,1) == 1 && ispk == 0     ispk = 1;    if Vall.chspk < 96        DATA.probelist = Vall.chspk;    else        DATA.probelist = AllV.GetProbeFromName(Vall.loadname);    end    newdata = 1;endif newdata > 0    DATA.duration = size(Vall.V,2) .* DATA.interval;    DATA.totalduration = DATA.duration;end