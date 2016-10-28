function h = MarkSelections(DATA, selected, calltag, cmenu, varargin)

oldh = getappdata(gcf,'Box3Handles');
for j =  1:length(oldh)
    if ishandle(oldh(j))
        delete(oldh(j));
    end
end
if size(selected,1) == length(DATA.exptlist)
    [a,b] = find(selected > 0);
    for j = 1:length(a)
        h(j) = PC.DrawBox(a(j),b(j),calltag , 'color','w','uicontextmenu',cmenu);
    end
    setappdata(gcf,'Box3Handles',h);
elseif iscell(selected)
    for j = 1:length(selected)
        h(j) = PC.DrawBox(selected{j}.exptid,selected{j}.probe,'','color','w');
    end        
    setappdata(gcf,'Box3Handles',h);
elseif size(selected,2) == 3
    for j = 1:size(selected,2)
    end
end
