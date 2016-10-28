function value = GetValue(DATA,type, varargin)
 %value = GetValue(DATA,type, varargin)
 % 'AllSPikes'
 % 'Clusters' [, 'withfits']
 % 'CellList'
 % 'currentcell'
 % 'shaperange'
 
 if isfigure(DATA)
     F = DATA;
     DATA = get(F,'UserData');
 end
 strargs = cell2cellstr(varargin);
 
value = []; 
 if strcmp(type,'AllSpikes')
     value = getappdata(DATA.toplevel,'AllSpikes');
     if isempty(value)
         value{length(DATA.exptlist),DATA.nprobes} = [];
     end
 elseif strcmp(type,'AutoList')
     value = getappdata(DATA.toplevel,'AutoCellList');
 elseif strcmp(type,'AutoClusters')
     value = getappdata(DATA.toplevel,'AutoClusters');
     if sum(strcmp('withfits',strargs))
         fits = getappdata(DATA.toplevel,'AutoClusterFits');
         value = AddFitFields(value,fits)
     end
 elseif strcmp(type,'Clusters')
     fits = [];
     if sum(strcmp('autolist',strargs))
         listtype = 'autolist';
     else
         listtype = PC.GetValue(DATA,'listtype');
     end
     if strcmp(listtype,'autolist')
         value = getappdata(DATA.toplevel,'AutoClusters');
         if sum(strcmp('withfits',strargs))
             fits = getappdata(DATA.toplevel,'AutoClusterFits');
         end
     else
         value = getappdata(DATA.toplevel,'Clusters');
         if sum(strcmp('withfits',strargs))
             fits = getappdata(DATA.toplevel,'ClusterFits');
         end
     end
     value = AddFitFields(value, fits);
 elseif strcmpi(type,'CellList')
     listtype = PC.GetValue(DATA,'listtype');
     if strcmp(listtype,'autolist') && isfield(DATA.autolist,'CellList')
         value = DATA.autolist.CellList;
     else
         value = DATA.CellList;
     end
 elseif strcmp(type,'currentcell')
     value = DATA.currentcell; %check for > 1 selected later...
 elseif strcmp(type,'Expts')
     value = getappdata(DATA.toplevel,'Expts');
 elseif strcmp(type,'listtype')
     value = 'normal';
       if isfield(DATA,'autolist') && isfield(DATA.autolist,'listtype')
           value = DATA.autolist.listtype;
       end
 elseif strcmp(type,'shaperange') %range of probes to use for shape correlation
     listtype = PC.GetValue(DATA,'listtype');
     if strcmp(listtype,'autolist') && isfield(DATA.autolist.corrproberange)
         value = DATA.autolist.corrproberange;
     elseif isfield(DATA,'ArrayConfig') && ~isempty(DATA.ArrayConfig)
         value = array.GetValue(DATA.ArrayConfig,'shaperange');
     else
         if DATA.probes == 24
             value = 5;
         else
             value = 2;
         end
         fpritnf('NO Array Config - using default proberange %d\n',value);
     end
 elseif strcmp(type,  'layoutfile')
     if sum(strcmp('last',strargs));
             value = strrep(DATA.defaultlayoutfile,'.layout','last.layout');         
     elseif isempty(DATA.layoutfile) || ~exist(DATA.layoutfile)
         value = DATA.defaultlayoutfile;
         if DATA.options.uselastlayout;
             value = strrep(DATA.defaultlayoutfile,'.layout','last.layout');
         end
     else
         value = DATA.layoutfile;
     end
 elseif strcmp(type, 'configfile')
     if sum(strcmp('last',strargs))
         value = strrep(DATA.defaultconfig,'.config','last.config');
     elseif isempty(DATA.configfile)
         value = DATA.defaultconfig;
     else
         value = DATA.configfile;
     end
 end

  function value = AddFitFields(value,fits)
      if ~isempty(fits) && iscell(value{1})
         for j = 1:length(fits)
             for k = 1:length(fits{j})
                 value{j}{k} = CopyFields(value{j}{k},fits{j}{k});
             end
         end
     elseif ~isempty(fits) && isstruct(value{1})
         for j = 1:length(fits)
             for k = 1:length(fits{j})
                 value{j} = CopySFields(value{j},k,fits{j}{k});
             end
         end         
     end
