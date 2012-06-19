function ConvertSpikeDir(path,type,varargin)

d = dir(path);
for j = 1:length(d)
    if regexp(d(j).name,'p[0-9]*t[0-9]*.mat')
        fprintf('%s\n',d(j).name);
        ConvertSpikeFile([path '/' d(j).name],type);
    end
end




function ConvertSpikeFile(name, type, varargin)

load(name);
if strcmp(type,'toint') && strcmp(class(Spikes.values),'double');
    maxv = max(Spikes.values(:));
    Spikes.values = int16(Spikes.values .* 16000/maxv);
    Spikes.maxv = maxv;
    save(name,'Spikes');
end