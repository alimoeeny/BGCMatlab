function FullV = ClearErrors(FullV,type, varargin)
%fullv.ClearErrors(FullV, type) removes old errors from FullV
%type 'flat' removes old errors warning of flat traces

if isdir(FullV)
    if nargin == 1
        type = 'all';
    end
    ClearDir(FullV, type, varargin);
    return;
end
    
if ~isfield(FullV,'errmsg') && ~isfield(FullV,'Errors')
    return;
end

if iscellstr(type)
    for j = 1:length(type)
        FullV = fullv.ClearErrors(FullV,type{j},varargin{:});
    end
    return;
end

if strcmp(type,'flat')
    fstr = 'Empty/Flat';
elseif strncmp(type,'missingprobe',12)
    fstr = 'Missing probes';
end
if isempty(fstr)
    fprintf('fullv.ClearErrors: Unrecognized error type %s',type);
    return;
end
oldid = find(CellToMat(strfind(FullV.errmsg,fstr)));
errfile = BuildFileName(FullV,'error');
if exist(errfile)
    X = load(errfile);
    X.Errors = FixErrors(X.Errors);
    errid = find([X.Errors.exptno] == FullV.exptno);
    eid = find(CellToMat(strfind({X.Errors(errid).s},fstr))); 
    if ~isempty(eid)
        X.Errors(errid(eid)) = [];
        try
            save(errfile,'-struct','X');
        catch ME
            CheckExceptions(ME);
        end
    end
else 
    X.Errors = [];
end

for j = 1:length(oldid)
    k = oldid(j);
    fprintf('Removing message from %s:%s\n',datestr(FullV.errdata(k).time),FullV.errmsg{k});
end
gid = setdiff(1:length(FullV.errmsg),oldid);
if isempty(gid)
    FullV = rmfields(FullV,{'errmsg' 'errdata'});
else
    FullV.errmsg = FullV.errmsg(gid);
    FullV.errdata = FullV.errdata(gid);
end

function ClearDir(name, type, varargin)

errfile = BuildFileName(name,'error');
X = load(errfile);
X.Errors = FixErrors(X.Errors);
for j = 1:length(X.Errors)
    if sum(strcmp(X.Errors(j).progname, {'fullv.Check' 'LoadFullV' 'BuildFullV' 'Spike2'}))
        clear(j) = 1;
    else
        clear(j) = 0;
    end
end
X.Errors = X.Errors(clear ==0);
save(errfile,'-struct','X');

