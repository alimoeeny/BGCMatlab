function F = FindFig(tag,strip)
%F = FindFig(tag) Find figure with tag
%%F = FindFig(tag, strip) finds figure where tag after removing strip
%%matches first arg

F = findobj(get(0,'children'),'flat','type','figure');
tags = get(F,'tag');
if nargin > 1 && ~isempty(strip)
    for j = 1:length(tags)
        tags{j} = strrep(tags{j},strip,'');
    end
end
id = find(strcmp(tag,tags));
if ~isempty(id)
    F = F(id);
else
    F = [];
end