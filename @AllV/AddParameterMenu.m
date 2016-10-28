function AddParameterMenu(F, callback, tag)%builds a menu that selects parameters for plotting spikes%Idea is that a parameter name added here will be converted into%a value by AllV.GetValues. No other work should be needed%    DATA = GetDataFromFig(F);    if ~isfield(DATA,'TemplateLabels')        cprintf('red','Figure %d Can''t get UserData\n',double(F));        return;    end    it = findobj(allchild(F),'flat','tag',tag,'type','uimenu');    delete(it);    hm = uimenu(F,'label',tag,'Tag',tag);    sm = uimenu(hm,'label','PC','tag','ParameterPCMenu');    for j = 1:12        s = sprintf('PC%d',j);        newtag = sprintf('%s%s',tag,s);        uimenu(sm,'Label',s,'callback',{callback, tag, s},'tag',newtag);    end    sm = uimenu(hm,'label','Template','tag','TemplateMenu');    for j = 1:length(DATA.TemplateLabels)        s = DATA.TemplateLabels{j};        newtag = sprintf('%sTmpl%d',tag,j);        uimenu(sm,'Label',s,'callback',{callback, tag, s},'tag',newtag);    end    for j = 1:2        newtag = sprintf('%sTemplatePC%d',tag,j);        s = sprintf('TemplatePC%d',j);        uimenu(sm,'Label',s,'callback',{callback, tag, s},'tag',newtag);    end        A = getappdata(DATA.toplevel,'AllTemplates');    for t = 1:length(A)        sm = uimenu(hm,'label',sprintf('Template %d',t),'tag',sprintf('TemplateMenu%d',t));        for j = 1:length(DATA.TemplateLabels)            s = sprintf('T%d.%s',t,DATA.TemplateLabels{j});            newtag = sprintf('%sTmpl%d',tag,j);            uimenu(sm,'Label',s,'callback',{callback, tag, s},'tag',newtag);        end    end            sm = uimenu(hm,'label','ADC','tag','TemplateADCMenu');    probes = unique(DATA.vpts(:,[1 3]));    probes = union(probes, DATA.chspk);    pts = unique(DATA.vpts(:,[2 4]));    for j = 1:size(probes,1)        for k = 1:size(pts);            s = sprintf('ADC %d:%d',probes(j),pts(k));            newtag = sprintf('%sADC%d:%d',tag,probes(j),pts(k));            uimenu(sm,'Label',s,'callback',{callback, tag, s},'tag',newtag);        end    end        strs = {};    for j = 1:length(DATA.chspk)        strs{j} = sprintf('energy%d',DATA.chspk(j));    end    strs = {strs{:} 'energysum' 'spkvar' 'ADC1' 'ADC2' 'ADC1dvdt' 'ADC2dvdt' 'ADC1dvdy1' 'ADC1dvdy2'  'ADC1csd' 'Centroid' 'CentroidSplit'};    for j = 1:length(DATA.chspk)        strs{end+1} = sprintf('dvrange%d',DATA.chspk(j));        strs{end+1} = sprintf('accel%d',DATA.chspk(j));    end    for j = 1:length(strs)        str = strs{j};        if strcmp(str,'ADC1')            str = sprintf('%d:%d',DATA.vpts(1,1),DATA.vpts(1,2));        elseif strcmp(str,'ADC2')            str = sprintf('%d:%d',DATA.vpts(2,1),DATA.vpts(2,2));        end                    newtag = sprintf('%s%s',tag,str);        uimenu(hm,'Label',str,'callback',{callback, tag, strs{j}},'tag',newtag);    end        