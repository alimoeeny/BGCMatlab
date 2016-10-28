function CheckLogs(name, varargin)
%AllV.CheckLogs(dir)  read RunAll.log and determines
%status of a parallel RunAllVpcs job
%AllV.CheckLogs(dir, 'history', d)
%       show only activity in the lst d days (default is 2)
%AllV.CheckLogs(dir, 'all') includes all dates

historymax = 2;  %go back two days by default
showneed = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'all',3)
        historymax = Inf;
    elseif strncmpi(varargin{j},'history',3)
        j = j+1;
        historymax = varargin{j};
    elseif strncmpi(varargin{j},'showneed',7)
        showneed = 1;
    end
    j = j+1;
end

if ~exist(name)
    fprintf('%s does not exist\n',name);
    return;
end
A = GetArrayConfig(name);
txt = scanlines([name '/RunAll.log']);
pid = find(CellToMat(regexp(txt,'Worker [0-9]* Processing')));
cid = find(strncmp('Completed',txt,9));
for j = length(cid):-1:1
    s = txt{cid(j)};
    finished(j) = GetExptNumber(s);
    finishdate(j) = datenum(s(end-20:end));
end
finished = finished';
started = [];
lines = {};
for j = length(pid):-1:1
    s = txt{pid(j)};
    e = GetExptNumber(s);
    d = datenum(s(end-20:end));
    id = find(cid > pid(j) & finished ==e);
    if length(started) < e || started(e) == 0
        started(e) = d;
        if ~isempty(id)
            lines{end+1} = sprintf('Expt%d %s to %s',e,datestr(d),datestr(finishdate(id(end))));
        else
            [newline, p] = ReadClusterLog(name,e,d);
            if ~isempty(newline)
                lines{end+1} = newline;
                if p(2) < length(A.id)
                    needid(length(lines)) = 1+p(2);
                end
            end
        end
        if ~isempty(lines)
            expts(length(lines)) = e;
            ts(length(lines)) = d;
            lineid(length(lines)) = pid(j);
        end
    end
end
cid = find(strncmp('Requested',txt,9));
requested = [];
for j = length(cid):-1:1
    s = txt{cid(j)};
    e = GetExptNumber(s);
    d = datenum(s(end-20:end));
    if e > length(requested) || requested(e) == 0
        requested(e) = d;
        if started(e) == 0  %not yet started by Runall. Check Log
            [newline, p] = ReadClusterLog(name,e,d);
            if ~isempty(newline)
                lines{end+1} = newline;
                if p(2) < length(A.id)
                    needid(length(lines)) = 1+p(2);
                else
                    needid(length(lines)) = 0;
                end
                expts(length(lines)) = e;
                ts(length(lines)) = d;
                reqlineid(length(lines)) = cid(j);
            end
        elseif d > started(e)  %this request not reached start yet
            fd = max(finishdate(find(finished ==e)));
            if d <  fd%did finish
                lines{end+1} = sprintf('Expt%d Requested %s Finished %s',e,datestr(d),datestr(fd));
                needid(length(lines)) = 0;
            else
                lines{end+1} = sprintf('Expt%d Requested %s not Started',e,datestr(d));
                needid(length(lines)) = 1;
            end
            expts(length(lines)) = e;
            ts(length(lines)) = d;
            reqlineid(length(lines)) = cid(j);
        else %d < started(e) means should have beed listed already
            
        end
    end
end

if showneed
    Expts = ReadExptDir(name);
    exptnos = GetExptNumber(Expts);
    missing = setdiff(exptnos,expts);
    d = now - historymax;
    for j = 1:length(missing)
        [newline, p] = ReadClusterLog(name,missing(j),d);
        if ~isempty(newline)
            lines{end+1} = newline;
            if p(2) < length(A.id)
                needid(length(lines)) = 1+p(2);
            end
            requested(missing(j)) = d;
        end
    end
end


if max(expts) > length(requested)
    requested(max(expts)) = 0;
end
[t, eid] = sort(requested);
for j = 1:length(eid)
    if now - t(j) < historymax
        id = find(expts == eid(j));
        if ~isempty(id)
            fprintf('%s\n',lines{id(end)});
        elseif showneed
%            fprintf('Need %s\n',lines{id(end)});
        end
    end
end

if showneed
    id = find(needid);
    fprintf('\nExpts Not COmpleted:\n');
    for j = 1:length(id)
        fprintf('%s\n',lines{id(j)});
    end
    Expts = ReadExptDir(name);
    exptnos = GetExptNumber(Expts);
    missing = setdiff(exptnos,expts);
    for j = 1:length(missing)
        fprintf('Expt%d not Requeted\n',missing(j));
    end
end


function [newline, p] = ReadClusterLog(name, e,d)
newline = [];
need = 0;
p = [0 0];
    ename = sprintf('%s/ClusterLogExpt%d.txt',name,e);
    etxt = scanlines(ename);
    x = dir(ename);
    estr = sprintf('E%dP',e);
    id = find(strncmp(estr,etxt,length(estr)));
    if ~isempty(id)
        p = sscanf(etxt{id(end)},'E%dP%d');
        if x.datenum < d
            newline = sprintf('Outdated Log (%s) for Expt%d working since %s',x.date,e,datestr(d));
        else
            newline = sprintf('Expt%d working since %s. Finished P%d',e,datestr(d),p(2));
        end
    end
