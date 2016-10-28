function CloseChildren(parent, varargin)
%CloseChildren(parent, closes figures associated with parent
% - idendifed via appdata'ParentFigure' in the Childrens
%CloseChildren(parent,'closeall') also closes figures that were separated by
%KeepFigure
%CloseChildren(parent,'tag', str) also closes figures that were separated by
closeall = 0;
j = 1;
args = {};
tag = '';
while j <= length(varargin)
    if strncmpi(varargin{j},'closeall',8)
        closeall = 1;
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        tag = varargin{j};
    end
    j = j+1;
end
%find hidden figure too
%figs = findobj('type','figure'); %does not find hidden figures
c = allchild(0);
types = get(c,'type');
fid  = find(strcmp('figure',types));
figs = c(fid);

if ischar(parent)
else
    figs = setdiff(figs,parent);
end
[~,id] = sort(double(figs));
figs = figs(id);
nc = 0;
childfigs = [];
for j = 1:length(figs)
    skip = 0;
    if ~isempty(tag) 
        x = get(figs(j),'tag');
        if isempty(strfind(x,tag))
            fprintf('Not Closing %s\n',x);
            skip = 1;
        end
    end
    x = getappdata(figs(j),'ParentFigure');
    if skip || isempty(x);
    elseif ischar(parent)
        if strcmp(x,parent)
            nc = nc+1;
            childfigs(nc) = figs{j};
            close(figs(j));
        end
    elseif x == parent || ismember(x,childfigs)
        if ~isappdata(figs(j),'KeepFigure') || closeall
            if figs(j) ~= parent
                nc = nc+1;
                childfigs(nc) = figs(j);
            end
        end
    end
end
if nc > 0
    close(childfigs);
end
