function str = CheckExceptions(res, label)
%CheckExceptions(res, label)  runs through a result file
% from AllVPcs, RunAllVpcs, RunAllGridFiles, and prints out 
% any results that were exceptions from try/catch
% also prints error stack if res is a single Exception struct
str = '';
nerr = 0;
showerr = 1;
outfid = 0;

if nargin < 2
    label = '';
elseif isnumeric(label) && ~isempty(fopen(label)) %an open pipe
    outfid = label;    
elseif strcmp(label,'-silent')
    showerr = 0;
    label = '';
end

if isstruct(res) 
    if isfield(res,'cls')
    CheckExceptions(res.cls);
    return;
    elseif isfield(res,'errs')
        CheckExceptions(res.errs);
        return;
    end
end
if iscell(res)
    nc = prod(size(res));
    for j = 1:nc
        newerr = CheckError(res{j},j);
        if newerr > 0
        nerr = nerr+newerr;
        end
    end
elseif isobject(res) %a single exception
    if showerr
    cprintf('red', 'Error %s\n',res.message);
    cprintf('red', 'From %s\n',res.identifier);
    end
    str = sprintf('%s\n%s\n',res.message, res.identifier);
    for j = 1:length(res.stack)
        s = sprintf('line %d %s in %s\n',res.stack(j).line,res.stack(j).name,res.stack(j).file);
        if showerr
            cprintf('blue','line %d %s in %s\n',res.stack(j).line,res.stack(j).name,res.stack(j).file);
        end
        str = [str s];
    end
    str = strrep(str,'\','/');
    if outfid
        fprintf(outfid,'Exception!!:%s',str);
    end
end
if nerr > 0 && ~isempty(label) %label given = use a popup
    warndlg(sprintf('%d Errors',nerr),[label  'Exceptions']);
end


function iserr = CheckError(R, id)
    iserr = 0;
    
    if iscell(R)
        for j = 1:length(R)
            CheckError(R{j},id);
        end
        return;
    end
if isfield(R,'errstate')
    s = [];
    line = [];
    if isfield(R,'filename')
        s = [s R.filename];
    end
    if isfield(R,'exptid')
        s = [s ' E' num2str(R.exptid)];
    end
    if isfield(R,'exptid')
        s = [s 'P' num2str(R.probe)];
    end
    if isfield(R,'datadir')
        s = [s R.datadir];
    end
    if id > 0
        s = [s num2str(id)];
    end
    if isempty(s)
        s = 'Unknown Filename';
    end
    [a,b] = fileparts(R.errstate.stack(1).file);
    line = sprintf(' :%s line %d',b,R.errstate.stack(1).line);
    mycprintf('errors','%s:%s%s\n',s,R.errstate.message,line);
    iserr = 1;
elseif isfield(R,'exc')
    if isfield(R,'name')
        fprintf('Error For %s:',R.name);
    else
        fprintf('Recorded Error'); 
    end
    CheckExceptions(R.exc);
end
