function LoadSelectedSpikes(DATA,eid,pid)    ts = now;    [a,b] = fileparts(DATA.name);    mnk = GetMonkeyName(DATA.name);    [c,d] = fileparts(a);    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');    for j = 1:length(eid);        e = eid(j);        p = pid(j);        xs = '';        if rem(DATA.exptid(e),1) > 0.001            xs = 'a';        end        name = [DATA.name '/Spikes/' mnk b '.p' num2str(p)  't' num2str(floor(DATA.exptid(e))) xs '.mat'];        fprintf('Reading %s at %s\n',name,datestr(now));        AllSpikes{e,p} = ReadSpikeFile(name);        AllSpikes{e,p}.probe = p;    end    fprintf('Spike Load took %.2f\n',mytoc(ts));    setappdata(DATA.toplevel,'AllSpikes',AllSpikes);    