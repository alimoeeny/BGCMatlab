function x = findobj(DATA, str, varargin)
%x = findobj(DATA, name) find gui elements for AllV
x = [];

if ~isfield(DATA,'toplevel')
    return;
end
if strcmp(str,'revert') && isfigure(DATA.toplevel)
    a = findobj(allchild(DATA.toplevel),'flat','tag','ClusterMenu');
    if ~isempty(a)
        b = findobj(allchild(a),'flat','tag','RevertCluster');
        x = findobj(allchild(b),'flat','label','Revert to');
    end
elseif strcmp(str,'revertlist') && isfigure(DATA.toplevel)
    a = findobj(allchild(DATA.toplevel),'flat','tag','ClusterMenu');
    if ~isempty(a)
        b = findobj(allchild(a),'flat','tag','RevertCluster');
        x = findobj(allchild(b),'flat','label','List Backups');
    end
elseif strcmp(str,'revertlist') && isfigure(DATA.toplevel)
elseif strcmp(str,'ChooseTrial') && isfigure(DATA.toplevel)
    F = findobj(allchild(DATA.toplevel),'flat','type','figure','Tag','Spikes');
    if ~isempty(F)
        x = findobj(allchild(F),'flat','tag','ChooseTrial');
    end
else
    x = findobj(allchild(DATA.toplevel),'flat','tag', str);
    if isempty(x)
        x = findobj(DATA.toplevel,'tag', str);
    end
end