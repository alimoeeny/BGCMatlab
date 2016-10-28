function errs = SaveErrors(FullV, varargin)
%take any errors in a FullV file and save them into a single
%file in each folder logging all FullV Errros

X.Errors= {};
silent = 0;
Errors.errs = {};
clearerrs = 0;
backup =1;
errs = {};
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'nerr') 
        Errors = varargin{j}; 
        if ~isfield(Errors,'errs')
            Errors.errs = {};
        end
    elseif strncmp(varargin{j},'clearerr',8)
        clearerrs = 1;
    elseif strncmp(varargin{j},'nobackup',6)
       backup = 0;
    end
    j = j+1;
end

if isfield(FullV,'errs') && ~isfield(FullV,'errmsg')
    FullV.errmsg = FullV.errs;
end

if isfield(Errors,'chstd')
    cpf = {'chstd' 'size' 'checkmains' 'chrange' 'checktime' 'nerr'};
    e = Errors.exptno;
    FullV.exptdetails(e).exptno = e;
    FullV.exptdetails(e).chstd = Errors.chstd;
    for k = 1:length(cpf)
        f = cpf{k};
        if isfield(Errors,f)
            FullV.exptdetails(e).(f) = Errors.(f);
        end
    end
end

%Load any existing Errors file
strs = {}; %list of previous errors
name = GetName(FullV,'folder');
if isdir(name)
    ename = [name '/Errors.mat'];
    if exist(ename)
        X = my.load(ename);
        if isfield(X,'Errors')
            if iscell(X.Errors)
                strs = CellToMat(X.Errors,'s');
            elseif isfield(X.Errors,'s')
                strs = {X.Errors.s};
            end
        end
        if ~isfield(X,'checktimes') && isfield(X,'exptdetails') && isfield(X.exptdetails,'checktime')
            for j = 1:length(X.exptdetails)
                X.checktimes(X.exptdetails(j).exptno) = X.exptdetails(j).checktime;
            end
        end
    end
else
    ename = '';
end

newerrs = 0;
if isfield(FullV,'exptdetails') 
    id = unique([FullV.exptdetails.exptno]);
    id = id(id>0);
    if ~isfield(X,'exptdetails') && ~isempty(id)
        X.exptdetails = FullV.exptdetails;
        newerrs = 2;
    elseif ~isempty(id)
        if isfield(X.exptdetails,'checktime') && length(X.checktimes) >= max(id)
            xid = find([FullV.exptdetails(id).checktime] > [X.checktimes(id)]);
            xid = id(xid);
        else
           xid = id;
        end
        if ~isempty(xid)
            X.exptdetails = CatStruct(X.exptdetails,FullV.exptdetails(xid),xid);
            X.checktimes(xid) = [FullV.exptdetails(xid).checktime];
            newerrs = 2;
        end
    end
end

if isfield(FullV,'checktimes') && isfield(X,'checktimes')
    if length(FullV.checktimes) > length(X.checktimes)
        xid = length(X.checktimes)+1:length(FullV.checktimes);
        X.checktimes(xid) = FullV.checktimes(xid);
        newerrs = 2;
    end
    xid = find(FullV.checktimes > X.checktimes(1:length(FullV.checktimes)));
    if ~isempty(xid)
        X.checktimes(xid) = FullV.checktimes(xid);
        newerrs = 2;
    end
elseif isfield(FullV,'checktimes')
    X.checktimes  = FullV.checktimes;
    newerrs = 2;
end

if (~isfield(FullV,'errmsg') || isempty(FullV.errmsg)) && isempty(Errors.errs) && ~clearerrs %nothing to do
    if isfield(FullV,'Errors')
        %nned to check
    elseif newerrs == 0
        return;
    end
end




if isfield(FullV,'olderrs')
    if ~isfield(X,'oldErrors')
        X.oldErrors = [];
    end
    rmid = [];
    for j = 1:length(FullV.olderrs)
        if iscellstr(FullV.olderrs(j).s)
            id = [];
            for k = 1:length(FullV.olderrs(j).s)
                a = find(strcmp(FullV.olderrs(j).s{k},strs));
                id = [id a];
            end
        else
            id = find(strcmp(FullV.olderrs(j).s,strs));
        end
        if ~isempty(id)
            eid = GetExptNumber(X.Errors(id));
            if isnan(eid) | eid == FullV.exptno | eid == 0
                fprintf('Removing Old Error %s\n',strs{id(1)});
                rmid = [rmid id];
            end
        end        
    end
    gid = setdiff(1:length(X.Errors),rmid);
    X.oldErrors = X.Errors(rmid);
    X.Errors = X.Errors(gid);
    strs = strs(gid);
    newerrs = 1;
end

if ~isfield(FullV,'errmsg')
    FullV.errmsg = {};
end
if ~isempty(FullV.errmsg) && length(FullV.errdata) < length(FullV.errmsg)
    FullV.errdata{length(FullV.errmsg)} = [];
end

for j = 1:length(FullV.errmsg)
    if isempty(X.Errors)
        X.Errors(1).s = FullV.errmsg{j};
        if iscell(FullV.errdata)
            X.Errors = CopySFields(X.Errors,1,FullV.errdata{j});
        else
            X.Errors = CopySFields(X.Errors,1,FullV.errdata(j));
        end
        newerrs = newerrs+1;
    elseif sum(strcmp(FullV.errmsg{j},strs)) == 0 %new error
        X.Errors(end+1).s = FullV.errmsg{j};
        if iscell(FullV.errdata)
            X.Errors = CopySFields(X.Errors,length(X.Errors),FullV.errdata{j});
        else
            X.Errors = CopySFields(X.Errors,length(X.Errors),FullV.errdata(j));
        end
        if ~silent
            fprintf('%s\n',FullV.errmsg{j});
        end
        newerrs = newerrs+1;
    end
end
if isfield(FullV,'Errors')
     for j = 1:length(FullV.Errors)
         if sum(strcmp(FullV.Errors(j).s,strs)) == 0 %new error
             X.Errors = CatStruct(X.Errors,FullV.Errors(j));
             newerrs = newerrs+1;
         end
     end
end

for j = 1:length(Errors.errs)
    if isempty(X.Errors)
        X.Errors(1) = Errors.errdata(j)
        X.Errors(1).s = Erros.errs{j};
    elseif sum(strcmp(Errors.errs{j},strs)) == 0 %new error
        X.Errors(end+1).s = Errors.errs{j};
        X.Errors = CopySFields(X.Errors,length(X.Errors),Errors.errdata(j));
    end
end



if newerrs
    good = [];
    for j = 1:length(X.Errors)
        if isempty(X.Errors(j).s)
            good(j) = 0;
        else
            good(j) = 1;
        end
    end
    X.Errors = X.Errors(good>0);
    if isfield(X,'Errors') && isfield(X.Errors,'s')
        [~,id,ia] = unique({X.Errors.s});
        [a,b] = Counts(ia);
        did = find(a>1);
        for k = 1:length(did)
            dup = find(ia == b(did(k)));
            id(did(k)) = dup(end);
        end
        X.Errors = X.Errors(id);
    end
    X.Errors = FixErrors(X.Errors);
    try
        if ~isempty(ename)
            if backup
                BackupFile(ename,'print');
            end
            save(ename,'-struct','X');
        else
            cprintf('red','!!!!\n!!!! Cant Save Errors for %s, becuse no folder %s\n',GetName(FullV),name);
        end
    catch ME
        CheckExceptions(ME);
    end
end

errs = X;