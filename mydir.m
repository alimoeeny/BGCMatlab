function [d, res] = mydir(path, pattern,varargin)
% d = mydir(path)
%works just like dir, but the name field contains full path.
%In order to avoid trouble with wildcards, if path is just hte
%directory name, have '/' at the end
%if path is a cell array of strings, calls mydir for each
% d = mydir(path, pattern) returns only entries that match regexp pattern
% d = mydir(path,'-dir') returns only directoriess
% d = mydir(path,'-dir','exclude',pattern) excludes names matching pattern
% (with regexp)
% d = mydir(path,'sortbyexpt') sorts result by expt number

exclude = '';
strargs = {};
args = {};
dironly = 0;
if nargin < 2
    pattern = [];
end
callback = [];
listtype = '';
callargs = {};
if isa(pattern,'function_handle')
    callback = pattern;
    pattern = [];
elseif iscell(pattern) && isa(pattern{1},'function_handle')
    callback = pattern{1};
    if length(pattern) > 1
        callargs = pattern(2:end);
    end
    pattern = [];
end
if strcmp(pattern,'-dir')
    dironly = 1;
    pattern = [];
    args = {args{:} '-dir'};
end
if strcmp(pattern,'-lrt')
    listtype = 'date';
    pattern = [];
end
if strcmp(pattern,'sortbyexpt')
    strargs{1} = pattern;
    pattern = [];
end
if ischar(pattern) && strcmp(pattern,'exclude') && ~isempty(varargin)
    exclude = varargin{1};
    pattern = [];
end

j = 1;
while j <= length(varargin)
    if isa(varargin{j},'function_handle')
        callback = varargin{j};
    elseif iscell(varargin{j}) && isa(varargin{j}{1},'function_handle')
        callback = varargin{j}{1};
        callargs = varargin{1}(2:end);
    elseif strncmpi(varargin{j},'exclude',5)
        j = j+1;
        exclude = varargin{j};
    elseif ischar(varargin{j})
        strargs = {strargs{:} varargin{j}};
    end
    j = j+1;
end

useregexp = 0;
if iscellstr(path)
    good = [];
    for j = 1:length(path)
        if ~isempty(pattern)
            d{j} = mydir([path{j} '/' pattern],varargin{:});
        else
            d{j} = mydir(path{j},args{:},varargin{:});
        end
        good(j) = ~isempty(d{j});
    end
    d = cat(1,d{find(good)});
else
    if isdir(path)
        root = path;
    else
        [root, name] = fileparts(path);
    end
    d = dir(path);
    if useregexp
    if isempty(d) %may have used regexp
        nx = 0;
        x  = dir(root);
        clear d;
        for j = 1:length(x)
            if regexp(x(j).name,name);
                nx = nx+1;
                x(j).filename = [root '/' x(j).name];
                d(nx) = x(j);
            end
        end
        return;
    end
    end
    if ~isempty(pattern)
        good = [];
        for j = 1:length(d)
            if iscellstr(pattern)
                good(j) = sum(CellToMat(regexp(d(j).name,pattern)));
            elseif regexp(d(j).name,pattern)
                good(j) = 1;
            end
        end
        d = d(find(good));
        clear good;
    end
    for j = 1:length(d)
        d(j).filename = d(j).name;
        if sum(strcmp(d(j).name,{'..' '.'}))
            good(j) = 0;
            else
            d(j).name = [root '/' d(j).name];
            good(j) = 1;
        end
        if dironly && ~isdir(d(j).name)
            good(j) = 0;
        end        
        if ~isempty(exclude) && ~isempty(regexp(d(j).name,exclude))
            good(j) = 0;
        end
    end
    if ~isempty(d)
    d = d(find(good));
    end
end

if sum(strcmp('sortbyexpt',strargs))
    eid = GetExptNumber({d.name});
    [a,b] = sort(eid);
    d = d(b);
elseif sum(strcmp('sortbydate',strargs)) || strcmp(listtype,'date')
    [a,b] = sort([d.datenum]);
    d = d(b);
end

if ~isempty(listtype)
    for j = 1:length(d)
        fprintf('%s %s %.2fkB\n',d(j).name,d(j).date,d(j).bytes./(1024));
    end
end

if isa(callback,'function_handle')
    for j = 1:length(d)
        fprintf('File %s\n',d(j).name);
        res{j} = feval(callback, d(j).name, callargs{:});
    end
end