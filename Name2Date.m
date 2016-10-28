function date = Name2Date(name, varargin)
%Name2Date(name) calculate date from name of file, assuming
%naming convention follows BackupFile
%
%see also BackupFile

%A special case with Pxx added to the end
if iscell(name)
    for j = 1:length(name)
        date(j) = Name2Date(name{j},varargin{:});
    end
    return;
elseif isfield(name,'loadname')
    date = Name2Date(name.loadname,varargin{:});
    return;
end

if ~ischar(name)
    date = NaN;
    return;
end

if regexp(name,'.*ClusterTimes([0-9,.]+)P[0-9]+.mat')
    dstr = regexprep(name,'.*ClusterTimes([0-9,.]+)P[0-9]+.mat','$1');
    dstr = [dstr(1:2) '-' dstr(3:4) '-' dstr(5:6) ' ' dstr(8:9) ':' dstr(10:11) '.' dstr(12:13)];
else    
    dstr = regexprep(name,'.*ClusterTimes([0-9,.]+).mat','$1');
    dstr = [dstr(1:2) '-' dstr(3:4) '-' dstr(5:6) ' ' dstr(8:9) ':' dstr(10:11) '.' dstr(12:13)];
end
try
    date = datenum(dstr);
    if nargout == 0
        fprintf('%s',datestr(date));
    end
catch ME
    CheckExceptions(ME);
end