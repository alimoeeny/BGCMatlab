function s = OldExLabel(DATA, j) s = ''; Clusters = getappdata(DATA.toplevel,'Clusters'); Expts = getappdata(DATA.toplevel,'Expts'); k = floor(Clusters{j}{1}.exptno); if DATA.show.exptno     s = [s num2str(k)]; end if DATA.show.exptname && ~isempty(Expts)     s = [s ':' Expts{j}.Header.expname]; end if DATA.show.ed && ~isempty(Expts)     s = [s ''  sprintf('%.2f', mean([Expts{j}.Trials.ed]))]; end