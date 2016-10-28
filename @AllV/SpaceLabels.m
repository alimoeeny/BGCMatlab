function str = SpaceLabels(DATA, C, space, varargin)
%str = SpaceLabels(DATA, space, varargin) return labels describing space


if ismember(space(1), [3 4])
    for j = 2:length(space)
        str{j-1} = DATA.TemplateLabels{space(j)};
    end
elseif space(1) == 1
    for j = 2:length(space)
        str{j-1} = sprintf('PC%d',space(j));
    end    
elseif space(1) == 2
    str = {};
    for j = 2:2:length(space)
        str{end+1} = sprintf('ADCP%dV%d',space(j),space(j+1));
    end    
elseif space(1) == AllV.USERSPACE
    str = {};
    for j = 1:length(C.Variables)
        str{end+1} = C.Variables{j};
    end    
elseif space(1) == 6
    str = {'AutoND' };
    if space(2) == 1
        str{2} = 'PC';
    end
else
    str{1} = 'Unknown';
end