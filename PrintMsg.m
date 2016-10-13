function s = PrintMsg(logfid, varargin)
%s = PrintMsg(logfid,format, args)   calls sprintf(format,args)
%then prints the string to the console, and into the file handle logfid
%(if logfid > 0)
%if logfid is a char and not empty, calls fopen(logfid,'a') then fclose.
%returns printed message

logname = '';
wasopen = 0;
s = sprintf(varargin{:});
fprintf('%s\n',s);
if ischar(logfid)
    if isempty(logfid)
        return;
    end
    logname = logfid;
    fid = fopen('all'); %check if this is already open
    for j = 1:length(fid)
        fname = fopen(fid(j));
        if strcmp(logname,fname)
            logfid = fid(j);
            wasopen = fid(j);
        end
    end
    if wasopen == 0
        logfid = fopen(logname,'a');
    end
else
    wasopen = 1;
end
if logfid > 0
    try
        fprintf(logfid,'%s\r\n',s);
        if wasopen == 0
            fclose(logfid);
        end
    catch
        fprintf('fid %d no longer valid\n',logfid);
    end
elseif ~isempty(logname)
    fprintf('Cant open log file %s\n',logname);
end
