function [sizes, f, xsizes, xf] = StructSize(S,varargin)
%find out which elements of a structure take the space

sizes = [];
xsizes = [];
xf = [];
verbose = 'default';
ts = now;
recurse = 1;
szf = {};

j = 1;
while j <=length(varargin)
    if strcmp(varargin{j},'silent')
        verbose = 'silent';
    elseif strcmp(varargin{j},'fields')
        j = j+1;
        szf = varargin{j};
    end
    j = j+1;
end
sz = whos('S');
b2m = (1024 .* 1024);
if ~isempty(szf)
    f = szf;
    for k = 1:length(szf)
        [sizes(k), xsizes(k)] = SumField(S,szf{k});
        fprintf('%.2fMb %s %d items\n',sizes(k)./b2m,szf{k},xsizes(k));
    end
    fprintf('Sum %.2f./%.2f\n',sum(sizes)./b2m,sz.bytes./b2m);
    return;
end

if iscell(S)
    allsize = [];
    allf = {};
    for j = 1:length(S)
        if isstruct(S{j}) || iscell(S{j})
            [sizes,f,xsizes,xf] = StructSize(S{j},'silent');
            for k = 1:length(f)
                if isfield(allsize,f{k})
                    allsize.(f{k}) = allsize.(f{k})+sizes(k);
                else
                    allsize.(f{k}) = sizes(k);
                end
            end
            for k = 1:length(xf)
                if isfield(allsize,xf{k})
                    allsize.(xf{k}) = allsize.(xf{k})+xsizes(k);
                else
                    allsize.(xf{k}) = xsizes(k);
                end
            end
        else
            X = S{j};
            a = whos('X');
            sizes(j) = a.bytes;
            if iscell(X)
                [sizes,f,xsizes,xf] = StructSize(S{j},'silent');
            elseif ischar(S{j})
                allsize.(genvarname(S{j})) = sizes(j);
            else
                allsize.(sprintf('arg%d',j)) = sizes(j);
            end
        end
        if mytoc(ts) > 10 && ~strcmp(verbose,'silent')
            fprintf('Finihsed Cell # %d/%d\n',j,length(S));
        end
    end
    if isempty(S)
        sizes(j) = 0;
        f = '';
    else
        f = fields(allsize);
        for j = 1:length(f)
            sizes(j) = allsize.(f{j});
        end
    end
elseif isstruct(S)
    f = fields(S);
    for j = length(f):-1:1
        X = rmfield(S,f{j});
        a = whos('X');
        sizes(j) = sz.bytes-a.bytes;
        if length(S) ==1 && iscell(S.(f{j})) && recurse
            [xsizes,xf] = StructSize(S.(f{j}),'silent');
            for k = 1:length(xf)
                xf{k} = [f{j} '_' xf{k}];
            end
        end
    end
else
    sizes(1) = sz.bytes;
    f{1} = 'var';
end
if ~strcmp(verbose,'silent')
    fprintf('Size(Mb)\tField\tTotal %.2fMb\n',sz.bytes./b2m);
    [ssize,sorted]  = sort(sizes);
    for j = 1:length(sorted)
        si = sorted(j);
        fprintf('%.3f\t%s (%.2f)\n',ssize(j)./b2m,f{si},sum(ssize(j:end))./b2m);
    end
end


function [sizes, ns] = SumField(S, f)

sizes = 0;
ns = 0;
if iscell(S)
    allsize = [];
    allf = {};
    for j = 1:length(S)
        if isstruct(S{j}) || iscell(S{j})
            [sz, n] = SumField(S{j},f);
            sizes = sizes+sz;
            ns = ns+n;
        end
    end
elseif isstruct(S)
    if isfield(S,f)
        sz = whos('S');
        X = rmfield(S,f);
        a = whos('X');
        sizes = sz.bytes-a.bytes;
        ns = ns+1;
    end
    sf = fields(S);
    for j = 1:length(sf)
        for k = 1:length(S)
            if isstruct(S(k).(sf{j})) || iscell(S(k).(sf{j}))
                [sz, n] = SumField(S(k).(sf{j}),f);
                sizes = sizes+sz;
                ns = ns+n;
            end
        end
    end
end