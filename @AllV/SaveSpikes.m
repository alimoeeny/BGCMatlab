function SaveSpikes(DATA, id, name)        if isempty(id)        return;    end    AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');    a = fileparts(name);    if ~exist(a,'dir')        mkdir(a);        fileattrib(a,'+w','g');    end    p = AllV.ProbeNumber(DATA);    Spikes.values = squeeze(AllVoltages(DATA.probe(1),:, id));    %thi sused to test for > 100, but that just messes up    %the file if ther eare a small number of spikes.    if size(Spikes.values,1) ~= length(DATA.spts)        cprintf('red','AllVoltages Size 2 is wrong\n');    end    %?check for values,1 = length(DATA.spts(    %changed in version 1.14    if size(Spikes.values,2) > 1        Spikes.values = Spikes.values';    end    Spikes.maxv = max(abs(Spikes.values(:)));    Spikes.maxint = double(intmax('int16')-5);    Spikes.VRange = [min(Spikes.values(:))  max(Spikes.values(:))];    Spikes.values = int16(Spikes.values .* Spikes.maxint./Spikes.maxv);    Spikes.times = reshape(DATA.t(id),length(id),1);    xy = DATA.xy{1}(id,1);    Spikes.codes = zeros(size(xy));    if isfield(DATA,'clst')        Spikes.codes = DATA.clst-1;    elseif DATA.cluster.sign < 0        Spikes.codes(xy < DATA.cluster.crit) = 1;    else        Spikes.codes(xy > DATA.cluster.crit) = 1;    end    Spikes.codes = reshape(Spikes.codes,length(Spikes.codes),1);    Spikes.Header.ctime = now;    Spikes.Header.version = DATA.version;    if isfield(DATA,'matfile');        Spikes.Header.matfile = DATA.matfile;    end    ts = now;    Spikes.Header.nspks = size(Spikes.values,1);    Spikes.Header.spts = DATA.spts;    PrintMsg(DATA.logname,'Saving %d Spikes to %s',length(id),name);    SpikeHeader = Spikes.Header;    save(name,'Spikes','SpikeHeader');    PrintMsg(DATA.logname,'Took %.2f\n',mytoc(ts));    if DATA.saveallspikes        xname = regexprep(name,'.p([0-9])*t','.p$1xt');%save channels needed to calculated dvdy also                [dspk, csdspk] = AllV.UseProbeList(DATA,'dvdy');        csdspk = union(dspk,csdspk);        chspk = union(DATA.chspk, csdspk);        chspk = setdiff(min(chspk):max(chspk), DATA.probe(1));        if ~strcmp(name,xname) && size(AllVoltages,1) > 1 && ~isempty(chspk)        Spikes.values = squeeze(AllVoltages(chspk,:, id));        Spikes.TriggerV = DATA.rV(id);        Spikes.maxv = max(abs(Spikes.values(:)));        Spikes.xVrange = [min(Spikes.values(:))  max(Spikes.values(:))];        Spikes.maxint = double(intmax('int16')-5);        Spikes.values = int16(Spikes.values .* Spikes.maxint./Spikes.maxv);        ts = now;        fprintf('Saving %d Other Probe Spikes to %s',length(id),xname);        Spikes.chspk = chspk;        save(xname,'Spikes');        fprintf('Took %.2f\n',mytoc(ts));        end    end    