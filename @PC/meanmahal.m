function d = meanmahal(DATA, id, a)Clusters = getappdata(DATA.toplevel,'Clusters');    expts = unique(DATA.allexpt(id));    for j = 1:length(expts)        d(j) = Clusters{expts(j)}{a}.mahal(1);    end    d = mean(d);