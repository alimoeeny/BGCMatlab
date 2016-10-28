function Figpos = GetFigPos(parent, varargin)
% GetFigPos(parent, varargin) find current locations of all child figures
% and then records these in appdata 'Figpos' in the parent
%if parent == 0, then stores all current figure locations in the root
%window appdata
%See also SetFigPos
Figpos = getappdata(parent,'Figpos');
f = findobj(get(0,'children'),'flat','type', 'figure');
apptag = '';
if isappdata(parent,'ApplicationTag')
    apptag = getappdata(parent,'ApplicationTag');
end

nc = 0;
childfigs = [];
[~,id] = sort(double(f));
f = f(id);

for j = 1:length(f)
    if double(parent) == 0
        p = 0;
    else
        p = getappdata(f(j), 'ParentFigure');
    end
    if p == parent | ismember(p,childfigs)
        nc = nc+1;
        tag = get(f(j),'tag');
        tag = regexprep(tag,'/','');
        if isfield(apptag,'extra')
            tag = strrep(tag,apptag.extra,'');
        end
        if isempty(tag)
            tag = spinrtf('Figure%d',double(f(j)));
        end  
        if isvarname(tag)
            Figpos.(tag) = get(f(j),'position');
        else
            Figpos.(genvarname(tag)) = get(f(j),'position');
        end
        childfigs(nc) = f(j);
    else
        ignorefig = f(j);
    end
end
if isempty(apptag)
    tag = genvarname(get(parent,'Tag'));
else
    tag = apptag.root;
end
if double(parent) > 0
    Figpos.(tag) = get(parent,'position');
end
setappdata(parent,'Figpos', Figpos);