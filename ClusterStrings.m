function dstr  = ClusterStings(C, strs, varargin)
if ischar(strs)
    str = strs;
    clear strs;
    strs{1} = str;
end
dstr = {};

for j = 1:length(strs)
    str = strs{j};
    if strcmpi(str,'mahal2D')
        dstr{j} = sprintf('%.1f',C.mahal(1));
    elseif strcmpi(str,'mahal1D')
        dstr{j} = sprintf('%.1f',C.mahal(4));
    elseif strcmpi(str,'mahalND')
        dstr{j} = sprintf('%.1f',C.mahal(1));
    elseif strcmpi(str,'fitdprime')
        if isfield(C,'fitdprime')
            dstr{j} = sprintf('%.1f',C.fitdprime(1));
        else
            dstr{j} = 'NA';
        end
    elseif strcmpi(str,'dropi')
        dstr{j} = sprintf('%.1f',C.dropi(3));
    elseif strcmpi(str,'date')
        if isfield(C,'savetime')
            dstr{j} = datestr(max(C.savetime),'mm/dd/yy');
        elseif isfield(C,'ctime')
            dstr{j} = ['*' datestr(max(C.ctime),'mm/dd/yy')];
        end
    elseif isfield(C,str)
        dstr{j} = num2str(C.(str));
    else
        dstr{j} = '';
    end
end