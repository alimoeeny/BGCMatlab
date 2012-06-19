function FullV = LoadFullV(name, varargin)
%FullV = LoadFullV(name)
%loads a FullV file from disk, and performs any necessary scaling
%

convert = 1;
FullV = [];
if ~exist(name,'file')
    fprintf('%s Does not exist\n',name);
    return;
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noconvert',5)
        convert = 0;
    end
    j = j+1;
end


ts = now;
load(name);
FullV.loadname = name;
method = 3;  %change manually for testing
FullV.initialloadtime = mytoc(ts);
d = whos('FullV');
FullV.readrate = d.bytes./FullV.initialloadtime;
FullV.size = d.bytes./(1024.*1024);
if convert && isinteger(FullV.V)  && isfield(FullV,'intscale')
    if method == 1  %try to save memory
        np = size(FullV.V,1);
        for j = 1:np
            NewV(np+1-j,:) = double(FullV.V(np+1-j,:)) .* FullV.intscale(1)/FullV.intscale(2);
            FullV.V = FullV.V(1:np-j,:);
        end
        FullV.V = NewV;
    elseif method == 2
        FullV.V = double(FullV.V) .*FullV.intscale(1)/FullV.intscale(2);
    elseif method == 3
        FullV.V = double(FullV.V);
        np = size(FullV.V,1);
        for j = 1:np
            FullV.V(j,:) = FullV.V(j,:) .* FullV.intscale(1)/FullV.intscale(2);
        end
    else
        NewV = double(FullV.V) .*FullV.intscale(1)/FullV.intscale(2);
        FullV.V = NewV;
        clear NewV;
    end
end
FullV.loadtime = mytoc(ts);

