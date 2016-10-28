function res = CallAllVPcs(DATA, eid, pid, varargin)  if isfigure(DATA)    DATA = GetDataFromFig(DATA);end   res = [];listtype = PC.GetValue(DATA,'listtype');cargs = {};xargs = {};  %args that need to be at the endj = 1;while j <= length(varargin)    if strncmpi(varargin{j},'forecfullv',6)        DATA.allvpcsmode = 'fullv'; %dont return and set - temporary    elseif strncmpi(varargin{j},'lastcut',6)        cargs = {cargs{:} 'applylast'};    elseif strncmpi(varargin{j},'autolist',3)        listtype = varargin{j};    elseif strncmpi(varargin{j},'noninteractive',6)        xargs = {xargs{:} 'noninteractive'};    elseif strncmpi(varargin{j},'refcut',3)        cargs = {cargs{:} 'refcut'};    elseif strncmpi(varargin{j},'setfit',6)        xargs = {xargs{:} varargin{j} varargin{j+1}};        j = j+1;    end    j = j+1;end    if DATA.usealltrials        args = {'usealltrials'};    else        args = {};    end    if strcmp(listtype,'autolist')        cargs = {cargs{:} 'useauto' 'autocutmode' 'ecker'};        CD = getappdata(DATA.toplevel,'AutoClusterDetails');        if length(CD) >= pid            if isempty(CD) || isempty(CD{eid})                C{1}.exptid = eid;                C{1}.probe = pid;                C{1}.exptno = DATA.exptlist(eid);                C = PC.AddFits(DATA,C);                CD{eid}{pid} = C{1};                CD{eid}{1}.exptno = DATA.exptlist(eid);            else                CD{eid} = PC.AddFits(DATA,CD{eid});            end            cargs = {cargs{:} 'ClusterDetails' CD{eid}};        end        cargs = {cargs{:} xargs{:}};    end    Expts = getappdata(DATA.toplevel,'Expts');    if length(Expts) >= eid && ~isempty(Expts{eid})        args = {args{:}, 'Expts', Expts};    end    ts = now;    AllFullV = getappdata(DATA.toplevel,'AllFullV');    if strncmp(DATA.allvpcsmode,'fromspikes',10)        name = sprintf('%s/Expt%dSpikes',DATA.name,DATA.exptid(eid));        args = {args{:} 'tchan' pid};        if strcmp(DATA.allvpcsmode,'fromspikesquick')            res = AllVPcs(name,args{:},'usecluster','nocheck', cargs{:});        else            res = AllVPcs(name,args{:},'reapply','nocheck', cargs{:});        end               elseif ~isempty(AllFullV)        exs = CellToMat(AllFullV,'exptno');        id = find(exs == DATA.exptid(eid));        args = {args{:} 'tchan' pid};        res = AllVPcs(AllFullV{id},args{:},DATA.allvpcsmode,'nocheck',cargs{:});    else                if strncmp(DATA.DataType,'Grid',4)            name = sprintf('%s/Expt%d.p%dFullV.mat',DATA.name,DATA.exptid(eid),pid);        else            name = sprintf('%s/Expt%dFullV.mat',DATA.name,DATA.exptid(eid));            args = {args{:} 'tchan' pid};        end        if exist(name)            res = AllVPcs(name,args{:},DATA.allvpcsmode,cargs{:});        end    end    clear CD; %takes time!! 0.5 sec.  Leave herer to avoid confusion    fprintf('Calling AllVPcs took %.1f sec\n',mytoc(ts));