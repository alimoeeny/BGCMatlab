function DATA = PlotMap(varargin)

%
%  PlotMap('monkeyname')
%  reads in penetration data for a monkey, and plots RF locations, and
%  the use grid
%
% PlotMap(fits)   plots data from rffits made by BuildRFFits
% PlotMap(fits, fixes) applies fixes (normally saved with the fits);
% PlotMap(dir)  if dir contains a file 'rffits.mat'  will load an plot it
%                      typically contains fits and fixes for an animal.
%
% Also reads penetration log summary. Some exceptions are handled by this
% Esp. if an array spans two areas.
%
% A line Boundary x Va Vb  says that above depth is area Va, below is Vb
% depth is caclulcated from ed + probe spacing
%
% See also BuildRFFits

name = 'RfMap';
% Set defaults before reading varargin
%
strings = [];
tag = 'RfMap'; %need to change these.
init = 0;
rffits = {};
MapDefs;

if length(varargin) & isnumeric(varargin{1})
    toplevel = varargin{1};
    varargin = varargin(2:end);
elseif length(varargin) & strncmpi(varargin{1},'new',3)
    toplevel = [];
    varargin = varargin(2:end);
else
    toplevel = findobj('Tag',tag);
end

if ~isempty(toplevel)
  if strncmpi(varargin{1},'store',5)
      DATA = varargin{2};
    set(toplevel,'UserData',varargin{2});
    return;
  else
    DATA = get(toplevel,'UserData');
  end
else
    DATA.monkeynames = {'duf' 'ruf' 'ica' 'lem' 'jbe'};
    DATA.useautomap = 0;
    DATA.tags.toplevel = tag;
    DATA.plot.area = 1;
    DATA.plot.labelpts = 1;
    DATA.plot.popuprf = 0;
    DATA.plot.type = 1;
    DATA.plot.fontsiz = 10;
    DATA.plot.maxage = 0;
    DATA.plot.minage = 0;
    DATA.verbose = 0;
    DATA.plot.auto = 0;
    DATA.plot.hemisphere = 2; %2 = Default Hemishpher for monkey
    DATA.plot.selecttype = 0;
    DATA.fittype = 'planar'; %seems best
    DATA.reloadpens = 0;
    prefsdir = '/b/group/matlab/preferences/PlotMap';
    DATA.layoutfile = [prefsdir '/' gethostname '.' GetUserName '.lastlayout.mat'];
    monkey = 'duf';
    init = 1;
end

j = 1;
if nargin 
    if ischar(varargin{1})
        if  strncmpi(varargin{1},'point',4)
            id = varargin{2};
            colors = mycolors;
            for j = 1:length(id)
                plotpen(DATA, DATA.px(id(j)), DATA.py(id(j)),'colors',colors);
                DATA.currentpen = [DATA.px(id(j)), DATA.py(id(j))];
            end
            type = gui.GetValue(DATA.top,'Tag','PenPlot');
            if strcmp(type,'PlotPen')
                ShowPenPlot(DATA,id(1));
            end
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'cellpoint',6)
            id = varargin{2};
            pe = DATA.map.pen(id,1);
            pname = sprintf('/b/data/ica/pens/pen%d.log',pe);
            GetFigure('OnePen','parent',DATA.toplevel);
            subplot(1,1,1);
            hold off;
            type = get(findobj('Tag','PenPlot'),'value') -1;
            if type == 4
                PlotOnePen(pname,'plot',type,'cmtype',[0 1 4]);
            elseif type > 0
                PlotOnePen(pname,'plot',type,'cmtype',1);
            end
        elseif  strncmpi(varargin{1},'cell',4)
            id = varargin{2};
            colors = mycolors;
            for j = 1:length(id)
                fprintf('%s %.2f %.2f %.1f\n',DATA.map.cellname{id(j)},DATA.map.rf(id(j),1),...
                    DATA.map.rf(id(j),2),DATA.map.depth(id(j)));
            end
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'excludepen',8)
            x = DATA.currentpen(1);
            y = DATA.currentpen(2);
            idx = find(DATA.map.pen(:,2) == x & DATA.map.pen(:,3) == y);
            DATA.excluded(idx) = 1;
            ReBuild(DATA);
        elseif  strncmpi(varargin{1},'missedpoint',8)
            id = find(DATA.map.pen(:,1) == varargin{2});
            fprintf('Missed on pen %d at %d,%d\n',varargin{2},DATA.map.pen(id,2), DATA.map.pen(id,3));
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'gridpoint',6)
            id = varargin{2};
            plotpen(DATA, DATA.px(id), DATA.py(id),'grid');
            DATA.currentpen = [DATA.px(id), DATA.py(id)];
            PlotMap(DATA.top,'store',DATA);
        elseif  strncmpi(varargin{1},'penpoint',6)
            id = find(DATA.map.pen(:,1) == varargin{2});
            pe = DATA.map.pen(id(1),1);
            fprintf('Pen %d\n',pe);
            pname = sprintf('/b/data/icarus/pens/pen%d.log',pe);
            GetFigure('OnePen','parent',DATA.toplevel);
            subplot(1,1,1);
            hold off;
            type = get(findobj('Tag','PenPlot'),'value') -1;
            if type == 4
                PlotOnePen(pname,'plot',type,'cmtype',[0 1 4]);
            elseif type > 0
                PlotOnePen(pname,'plot',type,'cmtype',1);
            end
        elseif  strncmpi(varargin{1},'depthpoint',4)
            ShowCell(DATA.map,varargin{2});
        elseif  strncmpi(varargin{1},'setmonkey',5)
            [DATA.map, DATA.monkey] = GetMap(DATA);
            DATA.selected = ones(size(DATA.map.area));
            DATA.excluded = zeros(size(DATA.map.area));
            DATA.marked = zeros(size(DATA.map.area));
            PlotMap(DATA.top,'store',DATA);
            RePlot(DATA);
        elseif  strncmpi(varargin{1},'savemap',7)
            SaveMap(DATA.map, DATA.map.monkeyname);
        elseif  strncmpi(varargin{1},'SetArea',5)
            DATA.plot.area = get(findobj('Tag','SetArea','Parent',toplevel),'value');
            set(toplevel,'UserData',DATA);
            PlotRFLoc(DATA);
        elseif  strncmpi(varargin{1},'getmap',5)
            return;
        elseif  strncmpi(varargin{1},'update',5)
            DATA = update(DATA);
            if DATA.plot.auto
                ReBuild(DATA);
            end
        elseif  strncmpi(varargin{1},'fittype',5)
            j = j+1;
            DATA.fittype = varargin{j};
            SetFitType(DATA);
        elseif  strncmpi(varargin{1},'plottype',5)
            DATA = update(DATA);
            ReBuild(DATA);          
        elseif  sum(strcmpi(varargin{1},DATA.monkeynames))
            monkey = varargin{j};
            init = 1;
            PlotMap(['/b/data/' monkey],varargin{:});
            return;
        elseif  strncmpi(varargin{1},'rufus',5)
            init = 1;
            monkey = 2;
        elseif  strncmpi(varargin{1},'icarus',5)
            init = 1;
            monkey = 3;
        elseif  strncmpi(varargin{1},'lem',3)
            init = 1;
            monkey = 4;
        elseif  strncmpi(varargin{1},'jbe',3)
            PlotMap('/b/data/jbe');
            return;
        elseif  strncmpi(varargin{1},'dae',3)
            init = 1;
            monkey = 5;
        elseif  strncmpi(varargin{1},'testmap',5)
            init = 1;
            monkey = 5;
        elseif strncmpi(varargin{1},'scroll',5)
            it = findobj('Tag','blockscroll','Parent',DATA.top);
            j = get(it,'value');
            scrolllist(j);
        elseif strncmpi(varargin{1},'checkblocks',5)
            [a, iset, values] = GetBlockList('Block');
            id = find(iset == 0);
            DATA.excluded(values(id)) = 1;
            PlotMap(DATA.top,'store',DATA);
            ReBuild(DATA);    
        elseif strncmpi(varargin{1},'clearblocks',5)
            DATA.excluded(find(DATA.excluded ==1)) = 0;
            PlotMap(DATA.top,'store',DATA);
            ReBuild(DATA);
        elseif strncmpi(varargin{1},'closeblocks',8)
            CloseTag('Blocks');
        elseif strncmpi(varargin{1},'choosepens',5)
            ChoosePens(DATA);
        elseif strncmpi(varargin{1},'markpen',5)
            x = DATA.currentpen(1);
            y = DATA.currentpen(2);
            idx = find(DATA.map.pen(:,2) == x & DATA.map.pen(:,3) == y);
            DATA.marked(idx) = 1;
            idx = find(DATA.xc == x & DATA.yc == y);
            DATA.penmarked(idx) = 1;
            PlotMap(DATA.top,'store',DATA);
        elseif strncmpi(varargin{1},'unmarkpen',5)
            x = DATA.currentpen(1);
            y = DATA.currentpen(2);
            idx = find(DATA.map.pen(:,2) == x & DATA.map.pen(:,3) == y);
            DATA.marked(idx) = 0;
        elseif strncmpi(varargin{1},'close',5)
            try
                SaveLayout(DATA,'nochoose');
            end
            CloseChildren(DATA.toplevel);
            close(DATA.top);
        elseif strncmpi(varargin{1},'plot',4)
            DATA = ReBuild(DATA);    
            RePlot(DATA);
        elseif strncmpi(varargin{j},'print',5)
            PrintPens(DATA,DATA.map);
        elseif isdir(varargin{1})
            mapfile = CheckNameBug([varargin{1} '/rffits.mat']);              
            if exist(mapfile);
                d = dir(mapfile);
                fprintf('Loading %s (%s)\n',mapfile,d.date);
                if length(varargin) > 1
                    args = varargin(2:end);
                else
                    args = {};
                end
                X = load(mapfile);
                X.loadname = mapfile;
                PlotMap(X,args{:});
                return;
            end
        else
            j = 1;
            init = 1;
            while(j < nargin)
                if(strncmpi(varargin{j},'name',3))
                    j = j+1;
                    name = varargin{j};
                end
                j = j+1;
            end
        end
    elseif iscell(varargin{1})
        DATA.map = rflist2map(DATA, varargin{j});
        rffits = varargin{j};
        if length(varargin) > j && isfield(varargin{j+1}{1},'fix')
            DATA.map = FixMap(DATA.map,varargin{j+1});            
        end
        DATA.monkey = DATA.map.monkey;
        DATA.map = ReadAllPens(DATA.map,DATA.map.monkeyname);
        monkey = DATA.monkey;
        
        init = 1;
    elseif isfield(varargin{1},'fits')
        X = varargin{1};
        args{1} = X.fits;
        if isfield(X,'fixes') && isfield(X,'extras')
            X.fixes = {X.fixes{:} X.extras{:}};
        end
        if isfield(X,'fixes')
            args = {args{:} X.fixes};
        end
        a = PlotMap(args{:},varargin{:});
        DATA = get(a.toplevel,'UserData');
        DATA.rffile = X.loadname;
        SetData(DATA);
        return;
    end
end

showgrid = 1;
j = 2;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'area') % a set of corrections
            DATA = FixMap(DATA,varargin{j});
        end
    elseif iscell(varargin{j})
        if isfield(varargin{j}{1},'fix')
            DATA.map = FixMap(DATA.map,varargin{j});
        end
    elseif strncmpi(varargin{j},'right',5) || strcmp(varargin{j},'R')
            DATA.plot.hemisphere = 1;
    elseif  strncmp(varargin{j},'nogrid',4)
        showgrid = 0;
    elseif  strncmp(varargin{j},'grid',4)
        showgrid = 1;
    elseif  strncmpi(varargin{j},'maxage',4)
        j = j+1;
        DATA.plot.maxage = varargin{j};
    elseif  strncmp(varargin{j},'MT',2)
        DATA.plot.area = MT;
    elseif  strncmp(varargin{j},'reload',4)
        DATA.reloadpens =1;
    end
    j = j+1;
end
if init & isempty(toplevel)
    scrsz = get(0,'Screensize');
    wsc = scrsz(3) /1280;  %scale everything as is 1280 wide
    if scrsz(3) > 2400
        wsc = 1.2;
    end
    wsize(1) = 380 * wsc;
    wsize(2) = 200 * wsc;
    
    
    if DATA.plot.hemisphere > 1
        DATA.plot.hemisphere = round(prctile(DATA.map.hemisphere,50));        
    end

    
    cntrl_box = figure('Position', [10 scrsz(4)-220*wsc 300*wsc 200*wsc],...
        'NumberTitle', 'off', 'Tag',tag,'Name',name);
    gui.SetDefaults(cntrl_box);
    DATA.top = cntrl_box;
    DATA.toplevel = cntrl_box;
    top = num2str(double(DATA.top));
    DATA.tags.fig = ['RFMapPlot' top];
    DATA.tags.figb = ['RFGridPlot' top];
    if( ~isempty(strings))
        lst = uicontrol(gcf, 'Style','listbox','String', strings,...
            'Callback', ' runcp(''setentry'')','Tag',listtag,...
            'Position',[10 10 wsize(1) wsize(2)]);
        
    end
    nr=7;
    HSPACE = 0.01;
    VSPACE = 0.01;
    cw = 0.02;
    nr = 6;
    ch = (1-nr.*VSPACE)./nr;
    bp(1) = HSPACE; bp(3) = 0.2; bp(2) = 1-ch; bp(4) = 0.9./nr;
    bp(3) = 0.23;
    row = 1;
    args = {'units' 'normalized'};
    xp = [0.01 1-row 0.1 1./nr];
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''next'');'],...
        'String', '>>', args{:},'Position', bp);
    xp(1) = xp(1) + xp(3);
    xp = [0.01 1-row 0.1 1./nr];
    bp(1) = bp(1) + bp(3) + HSPACE;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', [' PlotMap(' top ',''prev'');'],...
        'String', '<<',args{:}, 'Position', bp);
    
    bp(1) = bp(1) + bp(3) + HSPACE;
    uicontrol(gcf,'Style', 'pushbutton','String','Print','Callback',['PlotMap(' top ',''print'');'],args{:},'Position', bp);
    bp(1) = bp(1) + bp(3)+ HSPACE;
    uicontrol(gcf,'Style', 'pushbutton','String','Plot','Callback',['PlotMap(' top ',''plot'');'],args{:},'Position', bp);

    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 0.22;
    uicontrol(gcf,'Style', 'checkbox','String', 'Label Pts', 'Tag', 'LabelPts', 'Callback', [' PlotMap(' top ',''update'');'],...
        args{:},'Position', bp,'value',DATA.plot.labelpts);
    bp(1) = bp(1) + bp(3)+ HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'Pop RF', 'Tag', 'PopupRF', 'Callback', @UpdateChecks,...
        args{:},'Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'CheckBox','String','Verbose','Tag','Verbose',...
        'Value',DATA.verbose,'Callback',' PlotMap(DATA.top,''update'')',args{:},'Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'CheckBox','String','Auto','Tag','Auto',...
        'Value',DATA.plot.auto,'Callback',' PlotMap(DATA.top,''update'')',args{:},'Position', bp);
    
 %New row
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 0.3;
    bp(1) = HSPACE;
    uicontrol(gcf,'style','pop','string','Map|Grid|GridVar|GridDepth|Contour|Xpcolor|Ypcolor|dZdX|Arrows|MeanArrow|GMdepth|GMsmooth|Scatter|Fit|GridDate', ...
           'Callback', [' PlotMap(' top ',''plottype'');'], 'Tag','plottype',...
        args{:},'position',bp,'value',1); bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'text','String','Monkey',args{:},'Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    a = find(strcmp(monkey,DATA.monkeynames));
    uicontrol(gcf,'style','pop','string',DATA.monkeynames, 'Tag','RFMonkeyName','Callback',['PlotMap(' top ',''setmonkey'');'], args{:},'position',bp,'value',a);
  
    %New Row
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 0.1;
    uicontrol(gcf,'Style', 'text','String','Area',args{:},'Position', bp);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 0.15;
    uicontrol(gcf,'style','pop','string',areanames, 'Tag','SetArea','Callback',['PlotMap(' top ',''SetArea'');'], ...
        'value',DATA.plot.area,args{:},'position',bp);
    
    bp(1) = bp(1) + bp(3)+HSPACE;
    bp(3) = 0.2;
    uicontrol(gcf,'Style', 'text','String','HemiSphere',args{:},'Position', bp,'value',DATA.plot.hemisphere+1);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 0.2;
    uicontrol(gcf,'style','pop','string','L|R|Default','Tag','HemiSphere','Callback',['PlotMap(' top ',''update'');'],args{:}, 'position',bp,...
        'value',DATA.plot.hemisphere+1);

    bp(1) = bp(1) + bp(3)+HSPACE;
    bp(3) = 0.15;
    uicontrol(gcf,'Style', 'text','String','Fontsiz',args{:},'Position', bp);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 0.1;
    uicontrol(gcf,'style','pop','string','4|6|8|10|12', 'value',4,'Tag','FontSize','Callback',['PlotMap(' top ',''update'');'],args{:}, 'position',bp);
    

    
 %New Row
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 0.22;
    uicontrol(gcf,'Style', 'text','String','Max Age',args{:},'Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',DATA.plot.maxage),args{:},'Position', bp,'Tag','MaxAge',...
      'Callback', {@GuiMenu, 'maxage'},'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',0.0),args{:},'Position', bp,'Tag','MinAge', ...
	     'Callback', {@GuiMenu, 'minage'},'Backgroundcolor',[1 1 1]);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','None|Depths|RFs|3DRF|Comments|PlotPen', 'value',4,'Tag','PenPlot', args{:},'position',bp);

    
    
    
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    uicontrol(gcf,'Style', 'text','String','Select',args{:},'Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 0.2;
    uicontrol(gcf,'style','pop','string','All|MultContact', 'value',1,'Tag','SelectType','Callback',['PlotMap(' top ',''update'');'], args{:},'position',bp);
    
  hm = uimenu(gcf,'Label','File');
  top = DATA.top;
  uimenu(hm,'Label','Close','Callback',[' PlotMap(' num2str(double(DATA.top)) ',''close'');']);
  uimenu(hm,'Label','Reload','Callback',[' PlotMap(' num2str(double(DATA.top)) ',''setmonkey'');']);
  uimenu(hm,'Label','Choose','Callback',[' PlotMap(' double(top) ',''choosepens'');']);
  uimenu(hm,'Label','Exclude','Callback',[' PlotMap(' double(top) ',''excludepen'');']);
  uimenu(hm,'Label','Mark','Callback',[' PlotMap(' double(top) ',''markpen'');']);
  uimenu(hm,'Label','UnMark','Callback',[' PlotMap(' double(top) ',''unmarkpen'');']);
  uimenu(hm,'Label','Save Map','Callback',[' PlotMap(' double(top) ',''savemap'');']);
  hm = uimenu(gcf,'Label','Do');
  uimenu(hm,'Label','Refit','Callback',@MakeFits);
  uimenu(hm,'Label','List Pens','Callback',{@GuiMenu, 'list'});
  sm = uimenu(hm,'Label','Fit Type');
  uimenu(sm,'Label','Linear','tag','linear','callback',@SetFitType);
  uimenu(sm,'Label','Rotate','tag','rotate','callback',@SetFitType);
  uimenu(sm,'Label','Planar','tag','planar','callback',@SetFitType);
  uimenu(sm,'Label','Quadratic','tag','quadratic','callback',@SetFitType);
  uimenu(sm,'Label','Cubic','tag','cubic','callback',@SetFitType);
    set(gcf,'Menubar','none');
    if ~isfield(DATA,'map')
        [DATA.map, DATA.monkey] = GetMap(DATA);
    end
    if ~ischar(DATA.monkey)
        DATA.monkey = DATA.monkeynames{DATA.monkey};
    end
    DATA.selected = ones(size(DATA.map.area));
    DATA.excluded = zeros(size(DATA.map.area));
    DATA.marked = zeros(size(DATA.map.area));
    DATA.showgrid = showgrid;
    
    
    DATA = update(DATA);
    DATA.toplevel 
    DATA = ReBuild(DATA);
    set(DATA.toplevel,'UserData',DATA);
    if showgrid
        GetFigure(DATA.tags.figb,'parent',DATA.toplevel);
        hold off;
        PlotGrid(DATA);
    end
    ApplyLayout(DATA);
    if ~isempty(rffits)
        setappdata(DATA.toplevel,'rffits',rffits);
    end
end

function GuiMenu(src, event, type, varargin)

DATA = GetDataFromFig(src);
if strcmp(type,'list')
    MakeList(DATA);
elseif strcmp(type,'maxage')
    str = get(src,'string');
    DATA.plot.maxage  =  str2num(str);
    SetData(DATA);
elseif strcmp(type,'minage')
    str = get(src,'string');
    DATA.plot.minage  =  str2num(str);
    SetData(DATA);
end

function SetFitType(a,b);
if ishandle(a)
    DATA = GetDataFromFig(F);
    F = GetFigure(a);
    tag = get(a,'tag')
    DATA.fittype = tag;
else
    DATA = a;
end
DATA = PlotRFContour(DATA, 14,'diffs');
GetFigure(DATA.tags.figb,'parent',DATA.toplevel);
hold off;
PlotGrid(DATA);
SetData(DATA);

function UpdateChecks(a,b)
DATA = GetDataFromFig(a);
DATA = update(DATA);
SetData(DATA);

function DATA = MakeFits(a,b, varargin)

DATA = GetDataFromFig(a);
args = {};
if DATA.showarea == 1 %opercular V1
    args = {args{:} 'drange' 5};
end
DATA.fit = FitTopography(DATA.xc, DATA.yc, DATA.px, DATA.py,'rotate',args{:});
DATA.xyfit = FitTopography(DATA.xc, DATA.yc, DATA.px, DATA.py,'planar',args{:});
x(1,:) = [DATA.xyfit.xfit 0 0];
x(2,:) = [DATA.xyfit.yfit 0 0];
DATA.qfit = FitTopography(DATA.xc, DATA.yc, DATA.px, DATA.py,'quadratic','guess',x,args{:});
clear x;
x(1,:) = [DATA.qfit.xfit 0 0];
x(2,:) = [DATA.qfit.yfit 0 0];
DATA.cubicfit = FitTopography(DATA.xc, DATA.yc, DATA.px, DATA.py,'cubic','guess',x,args{:});
if strcmp(DATA.fittype,'planar')
    DATA.usefit = DATA.xyfit;
elseif strcmp(DATA.fittype,'quadratic')
    DATA.usefit = DATA.qfit;
elseif strcmp(DATA.fittype,'cubit')
    DATA.usefit = DATA.cubicfit;
else
    DATA.usefit = DATA.xyfit;
end
SetData(DATA);

function SetListAction(src, event)

F = GetFigure(src);
X = getappdata(F,'PlotOptions');
f = fields(X);
for j = 1:length(f)
    X.(f{j}) = 0;
end
tag = get(src,'tag');
onoff = get(src,'checked');
X.(tag) = 1;
SetMenuCheck(src, 'exclusive',X.(tag));
setappdata(F,'PlotOptions',X);

function MakeList(DATA)

[F, isnew] = GetFigure('PlotMapList','parent',DATA.toplevel);
if isnew
    gui.SetDefaults(F);
    set(F,'menubar','none');
    hm = uimenu(F,'Label','Action');
    X.plotall = 0;
    X.plotone = 1;
    X.plotrf = 0;
    X.plot3d = 0;
    X.showcommments = 0;
    uimenu(hm,'label','Plot All Matching pens','tag','plotall','callback',@SetListAction);
    uimenu(hm,'label','Plot This Pen','tag','plotone','callback',@SetListAction,'checked','on');
    uimenu(hm,'label','Plot This Pen 3D','tag','plot3d','callback',@SetListAction);
    uimenu(hm,'label','Plot Fits for This RF','tag','plotrf','callback',@SetListAction);
    uimenu(hm,'label','Show Commnnts','tag','showcomments','callback',@SetListAction);
    uimenu(hm,'label','Fix This Penetration','tag','fixpen','callback',@FixPen);
    bp = [0.02 0.02 0.96 0.96];
    lst = uicontrol(gcf,'Style', 'list','String','Max Age','Units','Normalized','Position', bp,'Tag','RFlist');
    setappdata(F,'PlotOptions',X);
    setappdata(F,'SelectedPen',1);
else
    lst = findobj(allchild(F),'flat','tag','RFlist');
end
map = DATA.map;

[idx, order] = sort(map.datenum)
nc = 0;
for j = 1:length(order)
    if DATA.selected(order(j))
        nc = nc+1;
        k = order(j);
        strs{nc} = sprintf('%s P%d (%.1f,%.1f) %s: %.1f,%.1f or%.0f',map.datestr{k},map.pen(k,1),...
            map.pen(k,2),map.pen(k,3),map.cellname{k},map.rf(k,1:2)-map.rf(k,9:10),map.rf(k,5));
        olist(nc) = order(j);
    else
        fprintf('Skipping P%d %.1f %.1f %s\n',map.pen(order(j),1),map.pen(order(j),2),map.pen(order(j),3),map.cellname{order(j)});
    end
end
set(lst,'String',strs,'callback',{@HitList},'UserData',olist);

function FixPen(src, event)

Parent = GetFigure(src);
pid = getappdata(Parent,'SelectedPen');
DATA = GetDataFromFig(Parent);
if ~isempty(pid)
    map = DATA.map;
    pe = map.pen(pid);
    P = DATA.map.pens{pe};
    rf = map.rf(pid,:);
end
[F, isnew] = GetFigure('Fix Penetration','parent',GetFigure(src));
if isnew
set(F,'menubar','none');
a = uicontrol(F,'style','Text','string','X','units','normalized','position',[0.01 0.01 0.1 0.5]);
b = gui.ui(F,'Edit',sprintf('%.1f',P.pos(1)),'left',a,'tag','PenXpos','height',0.5);
b = gui.ui(F,'Text','Y','left',b,'height',0.5);
b = gui.ui(F,'Edit',sprintf('%.1f',P.pos(2)),'left',b,'tag','PenYpos','height',0.5);
b = gui.ui(F,'TextEdit','Reason','','left',b,'width',0.4,'tag','Reason','height',0.5);
b = gui.ui(F,'pushbutton','Go','position',[0.85 0.01 0.1 0.15],'callback',@ApplyFix,'height',0.5);
b = gui.ui(F,'popupmenu',binoc.areanames,'up',a,'tag','VisualArea','height',0.5);
else
    a = gui.find(F,'PenXpos');
    set(a, 'string',sprintf('%.1f',P.pos(1)));
    b = gui.find(F,'PenYpos');
    set(b, 'string',sprintf('%.1f',P.pos(2)));
end
if isfield(P,'visualarea')
    a = unique(P.visualarea);
    id = find(strcmp(deblank(a{1}),binoc.areanames));
    if ~isempty(id)
        set(b,'value',id);
        b = gui.find(F,'VisualArea');
    end
end
setappdata(F,'SelectedPen',pid);

function ApplyFix(src, event)

Parent = GetFigure(src);
pid = getappdata(Parent,'SelectedPen');
DATA = GetDataFromFig(Parent);
areas = binoc.areanames;

if isfield(DATA,'rffile')
    xp = gui.Tag2Value('PenXpos',Parent);
    yp = gui.Tag2Value('PenYpos',Parent);
    X = my.load(DATA.rffile);
    pn = DATA.map.pen(pid);
    X.fixes{end+1}.pe = [pn xp yp];
    X.fixes{end}.fix = gui.Tag2Value('Reason',Parent);
    X.fixes{end}.name = DATA.map.cellname{pid};
    X.fixes{end}.area = gui.Tag2Value('VisualArea',Parent);
    X.fixes{end}.user = GetUserName;
    DATA.map = FixMap(DATA.map,X.fixes(end));
    DATA = ReBuild(DATA);
    DATA = RePlot(DATA);
    SetData(DATA);
    my.save(X);
end


function HitList(src, event);

DATA = GetDataFromFig(src);
F = GetFigure(src);
O = getappdata(F,'PlotOptions');
j = get(src,'value');
strs = get(src,'String');
s = strs{j};
order = get(src,'UserData');
j = order(j);
map = DATA.map;
pe = map.pen(j);
rf = map.rf(j,:);
fprintf('pen%d %s %.1f,%.1f\n',pe,map.cellname{j},rf(1:2));
DATA.plot.popuprf = 1;
x = GetFigure('OnePen','parent',DATA.toplevel);
hold off;
fits = getappdata(DATA.toplevel,'rffits');
if length(fits) >= map.fitid(j,1)
    fit = fits{map.fitid(j,1)};
    ifit = map.fitid(j,2);
else
    fit = [];
end
setappdata(F,'SelectedPen',j);

if O.plotone || O.plot3d
    P = DATA.map.pens{rf(6)};
    id = find(map.rf(:,6) == rf(6));
    for j = 1:length(id)
        P.rfs(j).depth = map.depth(id(j)).*1000;
        P.rfs(j).pos = map.rf(id(j),:);
        P.rfs(j).time = 0;
        P.rfs(j).name = sprintf('P%d',GetProbeNumber(map.cellname{id(j)}));
        if length(fits) >= map.fitid(id(j),1)
            fit = fits{map.fitid(id(j),1)};
            ifit = map.fitid(id(j),2);
            if isfield(fit,'fits')
                P.rfs(j).fits = fit.fits(fit.fitids{ifit});
            end
        end
    end
    if O.plot3d
        PlotOnePen(P,'rf3d');
    else
        PlotOnePen(P,'allcomments');
    end
elseif O.plotall
    plotpen(DATA, rf(7), rf(8),'colors',mycolors);
elseif O.showcomments
    s = [GetFilePath('data') '/' BuildFileName(DATA.map.cellname{j},'datadir')];
    s = regexprep(s,'P[0-9]+$','');
    s = [s '/'];
    PlotComments(s);
else
    fits = getappdata(DATA.toplevel,'rffits');
    if ~isempty(fits)
        GetFigure('RFPlots','parent',DATA.toplevel);
        hold off;
        PlotExptFit(fit.fits(fit.fitids{ifit}),'label');
    end
end



function area = Name2Area(s, varargin)
j = 1;
while j <= length(varargin)
    j = j+1;
end

areanums = [1 1 2 1 10 10 10 10 2]; %currently all unknown* cases are V1c
if strncmp(s,'Vc (',4)
    a = 5;
else
    a = find(strcmpi(s,{'V1' 'Vd' 'V2' 'unknown' 'Vc' 'unknown*' 'Vd*' 'V1*' 'V1Caclarine'}));
end
        if ~isempty(a)
            area = areanums(a);
        else
            area = 0;
        end

function map = rflist2map(DATA, R)
MapDefs;
areanums = [1 1 2 1 10 10 10 10]; %currently all unknown* cases are V1c
readpenlogs = 0;
hemisphere = 0;

for j = 1:length(R)
if isfield(R{j},'rf') || isfield(R{j},'rfs') || isfield(R{j},'proberf')
    good(j) = 1;
end
end
R = R(find(good));
nx = 0;

for j = 1:length(R)
    area = 1;
    if isfield(R{j},'name')
        name = R{j}.name;
    elseif isfield(R{j},'dirname')
        name = R{j}.dirname;
    end
    monk = GetMonkeyName(name);
    if isfield(R{j},'electrode') 
        if strcmp(R{j}.electrode,'uProbe')
            etype = MULTICONTACT;
        else
            etype = NORMAL;
        end
    else
        etype = NORMAL;
    end
    if isfield(R{j},'date')
        datestrs = datestr(R{j}.date);
        age = now-R{j}.date;
    else
        datestrs = '00/00/00';
        age = 0;
    end
    if isfield(R{j},'area')
        area = Name2Area(R{j}.area);
        map.area(j+nx) = area;
    end
    if isfield(R{j},'hemisphere')
    elseif strcmp(monk,'lem')
        hemisphere = 1;
    end
    if isfield(R{j},'proberf') && ~isempty(R{j}.proberf)  %has good fits
        for c = 1:size(R{j}.proberf,1)
            if c > 1
                nx = nx+1;
            end
            [mnk,mnkname,x,shortname] = GetMonkeyName(name);
            map.cellname{j+nx} = sprintf('%sP%d',regexprep(shortname,'/rffits.*',''),R{j}.proberf(c,12));
            map.rf(j+nx,:) = R{j}.proberf(c,1:10);
            map.depth(j+nx) = R{j}.proberf(c,11);
            map.probe(j+nx) = R{j}.proberf(c,12);
            map.area(j+nx) = area;
            map.hemisphere(j+nx) = hemisphere;
            monkey{j+nx} = monk;
            map.types(j+nx) = etype;
            map.age(j+nx) = age;
            map.datestr{j+nx} = datestrs;
            map.fitid(j+nx,1) = j;
            map.fitid(j+nx,2) = c;
        end
        etype = MULTICONTACT;
    elseif isfield(R{j},'rf') && isfield(R{j},'depth')
        map.rf(j+nx,:) = R{j}.rf;
        map.depth(j+nx) = R{j}.depth;
        map.cellname{j+nx} = regexprep(R{j}.name,'/rffits.*','');
        map.area(j+nx) = area;
        map.hemisphere(j+nx) = hemisphere;
        monkey{j+nx} = monk;
        map.types(j+nx) = etype;
        map.probe(j+nx) = 0;
        map.age(j+nx) = age;
        map.datestr{j+nx} = datestrs;
        map.fitid(j+nx,1) = j;
        map.fitid(j+nx,2) = 0;
    else
        monkey{j+nx} = monk;
        map.datestr{j+nx} = datestrs;
    end
end

map.datenum = datenum(map.datestr);
[a,b] = Counts(monkey);
[c,d] = max(a);
map.monkey = b{d};
map.monkeyname = b{d};
if readpenlogs
txt = scanlines(['/data/' map.monkey '/pens.err']);
for j = 1:length(txt)
    a = sscanf(txt{j},'%f');
    if length(a) > 2
        id = find(map.rf(:,6) == a(1));
        map.rf(id,7) = a(2);
        map.rf(id,8) = a(3);
    end
end
end
%.monkey is a char now 2015
%map.monkey = find(strcmp(map.monkeyname, DATA.monkeynames));
map.pen = map.rf(:,6:8);
if isfield(DATA,'toplevel') && isfigure(DATA.toplevel)
    setappdata(DATA.toplevel,'rffits',R);
end
    
function [map, monkey] = GetMap(DATA, varargin)
MapDefs;

global bgcfileprefix
monkey =  get(findobj('Tag','RFMonkeyName','Parent',DATA.top),'value');
strs = get(findobj('Tag','RFMonkeyName','Parent',DATA.top),'string');
monkey = strs{monkey};
nextfile = [];
if strncmpi(monkey,'duf',3)
    mapfile = '/bgc/bgc/anal/dufus/dufv1.fixtab';
    missfile = '/b/data/dufus/missed.pens';
    monkeyname = 'dufus';
elseif strncmpi(monkey,'ruf',3)
    mapfile = '/bgc//bgc/anal/rufus/rufv1.fixtab';
    missfile = '/b/data/rufus/missed.pens';
    monkeyname = 'rufus';
elseif strncmpi(monkey,'ica',3)
    mapfile = '/bgc//bgc/anal/icarus/icarus.fixtab';
    missfile = '/b/data/icarus/missed.pens';
    nextfile = '/b/data/icarus/penlist';
    monkeyname = 'icarus';
elseif strncmpi(monkey,'lem',3)
    mapfile = '/b/data/lem/lem.rftable';
    missfile = '/b/data/lem/missed.pens';
    nextfile = '/b/data/lem/penlist';
    monkeyname = 'lem';
elseif strncmpi(monkey,'dae',3)
    mapfile = '/b/data/dae/dae.fixtab';
    missfile = '/b/data/dae/missed.pens';
    nextfile = '/b/data/dea/penlist';
    monkeyname = 'dae';
    
elseif strncmpi(monkey,'test',4)
    mapfile = '/bgc//bgc/anal/maps/test.tab';
    missfile = '/b/data/lem/missed.pens';
else
    mapfile = sprintf('/b/bgc/anal/%s/%s.fixtab',monkey,monkey);
    missfile = sprintf('/b/bgc/anal/%s/missed.pens',monkey);
    nextfile = sprintf('/b/data/%s/penlist',monkey);
end
automapfile = strrep(mapfile,'.rftable','.fixtab');
if DATA.useautomap || (~exist(mapfile) && exist(automapfile))
    mapfile = automapfile;
end
DATA.monkey = monkey;
if exist('bgcfileprefix','var')
    %    mapfile = [bgcfileprefix mapfile];
end

map = loadmap(mapfile);
map.monkeyname = monkey;
if isfield(DATA,'map') & isfield(DATA.map, 'pens')
    map.pens = DATA.map.pens;
    map.missed = DATA.map.missed;
else
    args = {};
    if DATA.reloadpens
        args = {args{:} 'reload'};
    end
    map = ReadAllPens(map,map.monkeyname,args);
end

if monkey == RUFUS
    id = find(map.datenum < datenum('9/1/2006'));
    map.hemisphere(id) = 0;
    id = find(map.datenum > datenum('9/1/2006'));
    map.hemisphere(id) = 1;
end
if exist(missfile) & isnumeric(map.missed) & sum(map.missed) == 0
    [a,b] = textread(missfile,'%d %s');
    v = strmatch('V2',b);
    if ~isempty(v)
        missed{2} = a(v);
    end
    if isempty(a)
        missed{1} = [];
        missed{2} = [];
    end
        map.missed = missed;
end
if exist(nextfile)
    fid = fopen(nextfile,'r');
    a = textscan(fid,'%s','delimiter','\n');
    txt = a{1};
    for j = 1:length(txt)
        crd(j) = textscan(txt{j},'%n');
        nf(j) = length(crd{j});
        if nf(j) <= 2 | crd{j}(3) == 0
            unused(j) = 1;
        else
            unused(j) = 0;
        end
    end
    id = find(unused)
    if length(id)
        for j = 1:length(id)
        map.nextpens(1:length(crd{id(j)}),j) = crd{id(j)};
        end
    end
end

function DATA = update(DATA, varargin)
  DATA.plot.labelpts = get(findobj('Tag','LabelPts','Parent',DATA.top),'value');
  DATA.plot.popuprf = get(findobj('Tag','PopupRF','Parent',DATA.top),'value');
  DATA.plot.type = get(findobj('Tag','plottype','Parent',DATA.top),'value');
  DATA.plot.area = get(findobj('Tag','SetArea','Parent',DATA.top),'value');
  DATA.plot.selecttype = get(findobj('Tag','SelectType','Parent',DATA.top),'value');
  DATA.plot.hemisphere = get(findobj('Tag','HemiSphere','Parent',DATA.top),'value')-1;
  str = get(findobj('Tag','FontSize','Parent',DATA.top),'string');
  j = get(findobj('Tag','FontSize','Parent',DATA.top),'value');
  DATA.plot.fontsiz = str2num(str(j,:));
  DATA.verbose = get(findobj('Tag','Verbose','Parent',DATA.top),'value');
  DATA.plot.auto = get(findobj('Tag','Auto','Parent',DATA.top),'value');

  
  DATA.plot.maxage = str2num(get(findobj('Tag','MaxAge','Parent',DATA.top),'string'));
  DATA.plot.minage = str2num(get(findobj('Tag','MinAge','Parent',DATA.top),'string'));

  PlotMap(DATA.top,'store',DATA);

  
function CloseTag(tag)
it = findobj('Tag',tag);
if ~isempty(it)
    close(it);
end


function DATA = ReBuild(DATA);
    DATA = PlotRFLoc(DATA);
    
function DATA = RePlot(DATA)

MapDefs;
if DATA.plot.area == COMPAREV1V2
    gridlines = 0;
else
    gridlines = 1;
end
    
GetFigure(DATA.tags.fig,'parent',DATA.toplevel);
hold off;
if DATA.plot.type == 2
    GetFigure(DATA.tags.figb,'parent',DATA.toplevel);
    hold off;
    PlotGrid(DATA);
    return;
elseif DATA.plot.type == 3
    PlotGridVar(DATA);
    return;
elseif DATA.plot.type == 9
    PlotGridArrows(DATA);
    return;
elseif DATA.plot.type == 10
    PlotGridArrows(DATA,'mean');
    return;
elseif DATA.plot.type == 11
    PlotGMDepths(DATA);
    return;
elseif DATA.plot.type == 12
    PlotGMDepths(DATA, 'smooth', 0.3);
    return;
elseif DATA.plot.type == 4
    PlotGridDepths(DATA);
    return;
elseif ismember(DATA.plot.type, [5 6 7 13 14]) %contour/pcolor plot
    DATA = PlotRFContour(DATA, DATA.plot.type);
elseif DATA.plot.type == 8
    PlotZX(DATA);
elseif DATA.plot.type == 13
    PlotRFScatter(DATA)
elseif DATA.plot.type == 15
    PlotGridImage(DATA,DATA.plot.type);
    return;
else
    PlotRFMap(DATA, gridlines);
    offset = range(DATA.xc)/60;
    if isfield(DATA.map,'nextpens') & length(DATA.pid) > 10
        recent = (length(DATA.pid)-10):length(DATA.pid);
        id = find(~ismember(DATA.map.nextpens(1,:)+i*DATA.map.nextpens(2,:),DATA.ppx(recent)+i*DATA.ppy(recent)));
        for j = id;
        x = PredictRF(DATA, DATA.map.nextpens(:,j));
        plot(x(1),x(2),'go');
        text(x(1)+offset,x(2),sprintf('%d,%d',DATA.map.nextpens(1,j),DATA.map.nextpens(2,j)),'color','g');
        fprintf('Pen at %.1f,%.1f RFs should be %.1f,%.1f\n',DATA.map.nextpens(1,j),DATA.map.nextpens(2,j),x(1),x(2));
        end
    end
    if DATA.showgrid
        GetFigure(DATA.tags.figb,'parent',DATA.toplevel);
        hold off;
        PlotGrid(DATA);
    end

end
 
function px = PredictRF(DATA, x)
sd = 1;
C = (x(1) +i * x(2)) - (DATA.px + i * DATA.py);
w = exp(-(abs(C).^2)/2.*sd.^2);
px(1) = sum(w.*DATA.xc)./sum(w);
px(2) = sum(w.*DATA.yc)./sum(w);

function PlotZX(DATA)
plotdz = 6;
usefits = 1;
top = num2str(DATA.top);
pens = unique(DATA.map.pen(:,1));
dzs = [];
dxs = [];
dys = [];
deccs = [];
dangles = [];
ids = [];
for j = 1:length(pens)
    id = find(DATA.map.pen(:,1) == pens(j));
    dz = DATA.map.depth(id);
    [mz, mzi] = min(dz);
    dz = dz -mz;
    if usefits
        dx = DATA.map.rf(id,6)- DATA.map.rf(id(mzi),6);
        dy = DATA.map.rf(id,7)- DATA.map.rf(id(mzi),7);
        ecc = abs(DATA.map.rf(id,6)+ i * DATA.map.rf(id,7));
        angles = atan2(DATA.map.rf(id,6),DATA.map.rf(id,7));
    else
        dx = DATA.map.rf(id,1)- DATA.map.rf(id(mzi),1);
        dy = DATA.map.rf(id,2)- DATA.map.rf(id(mzi),2);
        ecc = abs(DATA.map.rf(id,1)+ i * DATA.map.rf(id,2));
        angles = atan2(DATA.map.rf(id,2),DATA.map.rf(id,1));
    end
    decc = ecc - ecc(mzi);
    dangle = angles - angles(mzi);
    for k = 1:length(id)
    if plotdz == 1
    h = plot(dz(k),dangle(k),'o');
    hold on;
    elseif plotdz == 6
        h = plot(dz(k),decc(k),'o');
        hold on;
    elseif plotdz == 5
        h = plot(dz(k),dx(k)-dy(k),'o');
        hold on;
    else
    h = plot(dz(k),dx(k),'o');
    hold on;
    h = plot(dz(k),dy(k),'ro');
    end
    set(h,'buttondownfcn',['PlotMap(' top ',''cell'',' num2str(id(k)) ',1);']);
    end
    if length(dz) > 1
        dangles = [dangles dangle(2:end)'];
        deccs = [deccs decc(2:end)'];
        ids = [ids id(2:end)'];
        dzs = [dzs dz(2:end)];
        dxs = [dxs dx(2:end)'];
        dys = [dys dy(2:end)'];
    end
end
[x,y] = meshgrid([min(dxs):0.01:max(dxs)],[min(dys):0.01:max(dys)]);
Z = Interpf(dxs, dys, dzs,x,y,1,0.01);
hold off;
if plotdz == 3
imagesc(x(1,:),y(:,1),Z);
elseif plotdz == 4
    plot3(dxs,dys,dzs,'o');
end
    

function ShowPenPlot(DATA,id)
pe = DATA.id(id);
pname = sprintf('/b/data/%s/pens/pen%d.log',DATA.monkey,pe);
GetFigure('OnePen','parent',DATA.toplevel);
subplot(1,1,1);
hold off;
type = gui.GetValue(DATA.toplevel,'Tag','PenPlot');
if strcmp(type,'Comments')
    PlotOnePen(pname,'plot',type,'cmtype',[0 1 4]);
elseif ~strcmp(type,'None')
    PlotOnePen(pname,'plot',type,'cmtype',1);
end


function DATA = PlotRFLoc(DATA, varargin)

MapDefs;
showall = 0;
holdon = 0;
map = DATA.map;
GetFigure(DATA.tags.fig,'parent',DATA.toplevel);
area = DATA.plot.area;

j = 1;
while(j < nargin)
    if strncmpi(varargin{j},'hold',4)
        holdon = 1;
    elseif strncmpi(varargin{j},'all',3)
        showall = 1;
    elseif strncmpi(varargin{j},'area',3)
        area = varargin{j+1};
        j = j+1;
    end
    j = j+1;
end

gridlines = 1;
if area == COMPAREV1V2
    gridlines = 0;
    idx = find(map.area ~=1);
    idx = find(map.area == V2);
    idxall = [];
    for j = idx
        v1 = find(map.area' == V1 & map.pen(:,2) == map.pen(j,2) & map.pen(:,3) == map.pen(j,3));
        idxall = [idxall v1'];
    end
    DATA.selected = zeros(size(map.area));
    DATA.selected(unique(idxall)) = 1;
    DATA.selected(idx) = 1;
elseif area ~= ALL_AREAS
    idx = find(map.area == area & map.hemisphere == DATA.plot.hemisphere);
    DATA.selected = zeros(size(map.area));
    DATA.selected(idx) = 1;
    if area <= length(map.missed)
        if iscell(map.missed)
            DATA.missed = find(map.missed{area} == area);
        else
            DATA.missed = find(map.missed == area);
        end
    else
        DATA.missed = [];
    end
else
    idx = 1:length(map.area);
    DATA.selected = ones(size(map.area));
    DATA.missed = [map.missed{1} map.missed{2}];
end

if DATA.monkey == DUFUS
    idx = find(map.pen(:,1) < 2);
else
    idx = find(map.pen(:,1) < 1);
end

    
DATA.selected(idx) = 0;
idx = find(DATA.excluded ==1);
DATA.selected(idx) = 0;

if DATA.plot.maxage > 0
    id = find(map.age > DATA.plot.maxage);
    DATA.selected(id) = 0;
end
if DATA.plot.minage > 0
    id = find(map.age < DATA.plot.minage);
    DATA.selected(id) = 0;
end

if DATA.plot.selecttype ==MULTICONTACT
    id = find(map.types ~= MULTICONTACT)
    DATA.selected(id) = 0;
end
    
PlotMap(DATA.top,'store',DATA);
DATA.showarea = area;

if DATA.plot.type == 4
    PlotGridDepth(DATA);
    return;
end

cellsonly = 1;

idx = find(DATA.selected > 0);
selected = idx;
if isempty(idx)
    fprintf('No data for Area %d\n',area);
    return;
end
tic;

for j = idx;
    xpos(j) = map.rf(j,1)-map.rf(j,9);
    ypos(j) = map.rf(j,2)-map.rf(j,10);
    if map.pen(j,2) > -99 & DATA.selected(j)
        all.px(j) = map.pen(j,2);
        all.py(j) = map.pen(j,3);
        all.area(j) = map.area(j);
        all.id(j) = map.pen(j,1);
        if regexp(map.cellname{j},'M[0-9][0-9][0-9]')
            all.ptype(j) = 1;
        elseif regexp(map.cellname{j},'S[0-9][0-9][0-9]')
            all.ptype(j) = 1;
        else
            all.ptype(j) = 0;
        end
    else
        all.px(j) = NaN;
        all.py(j) = NaN;
        all.area(j) = NaN;
        all.id(j) = NaN;
    end
    if strcmp(map.cellname{j},'Nocell')
        all.cell(j) = 0;
    else
        all.cell(j) = 1;
    end
end


toc;
idx = find(DATA.selected == 0);
all.px(idx) = NaN;
all.py(idx) = NaN;
all.cell(idx) = NaN;
all.id(idx) = NaN;

if ~holdon
    hold off;
end

if showall
    plot(xpos,ypos,'o');
    hold on;
end
np = 1;



%Collect data together by penetration location
DATA.xc = [];
DATA.yc = [];
DATA.px = [];
DATA.py = [];
DATA.rvar = [];
DATA.id = [];
DATA.ptype = []
%ids = DATA.id;
xvals = all.px(find(~isnan(all.px)));
yvals = all.py(find(~isnan(all.py)));
tic;
allzs = [];
for x= unique(xvals)
    rowstart = np;
    for y= unique(yvals);
        idx = find(all.px == x & all.py == y & all.cell > 0);
        if isempty(idx) && ~cellsonly
            idx = find(all.px == x & all.py == y);
        end
        if(length(idx) > 0)
            if length(idx) > 3 %remove outliers
                xzs = abs(zscore(xpos(idx)));
                yzs = abs(zscore(ypos(idx)));
                badid = find(xzs > 3 | yzs > 3);
                allzs = [allzs xzs(:)'];
                if ~isempty(badid)
                    for k = 1:length(badid)
                        a = idx(badid(k));
                        fprintf('Excluding %.2f,%.2f(%.2f) %s P%d\n',xpos(a),ypos(a),xzs(badid(k)),DATA.map.cellname{a},DATA.map.probe(a));
                    end
                    idx(badid) = [];
                end
            end
            if area == COMPAREV1V2
                v1 = find(all.area(idx) == 1);
                v2 = find(all.area(idx) ~= 1);
                DATA.xc(np) = mean(xpos(idx(v1)));
                DATA.yc(np) = mean(ypos(idx(v1)));
                DATA.axc(np) = mean(xpos(idx(v2)));
                DATA.ayc(np) = mean(ypos(idx(v2)));
            else
                DATA.xc(np) = mean(xpos(idx));
                DATA.yc(np) = mean(ypos(idx));
            end
                rad = sqrt((xpos(idx) - DATA.xc(np)) .^2 + (ypos(idx) - DATA.yc(np)).^2);
            if length(rad) > 1
                DATA.rvar(np) = std(rad);
            else
                DATA.rvar(np) = 0;
            end
%            DATA.id(np) = ids(idx);
            if sum(DATA.marked(idx)) > 0
                DATA.penmarked(np) = 1;
            else
                DATA.penmarked(np) = 0;
            end      
            DATA.id(np) = max(all.id(idx));
            DATA.allid{np} = unique(all.id(idx));
                
            DATA.px(np) = x;
            DATA.py(np) = y;
            if sum(map.types(idx) == MULTICONTACT) % some multicontact probes 
                DATA.ptype(np) = MULTICONTACT;
            else
                DATA.ptype(np) = 0;
            end
            np = np+1;
            if(DATA.verbose)
                fprintf('%.0f %.0f %d\n',x,y,np);
            end
        end
    end
%    plot(DATA.xc(rowstart:np-1),DATA.yc(rowstart:np-1),':');
end
np = 0;
for pid= unique(all.id)
    rowstart = np;
    idx = find(all.id == pid);
        if(length(idx) > 0) & pid> 0 & sum(all.cell(idx)) > 0
            np = np+1;
                DATA.pxc(np) = mean(xpos(idx));
                DATA.pyc(np) = mean(ypos(idx));
                DATA.pid(np) = pid;
            DATA.ppx(np) = mean(all.px(idx));
            DATA.ppy(np) = mean(all.py(idx));
            pxs = unique(all.px(idx));
            pys = unique(all.py(idx));
            if length(pxs) >1 | length(pys) > 1
                fprintf('Pen %d 2 Co-ordinates: px ',pid);
                fprintf('%.1f,',pxs);
                fprintf(' py: ');
                fprintf('%.1f,',pys);
                fprintf('\n');
            end
            DATA.ppid(np) = find(DATA.px == all.px(idx(1)) & DATA.py == all.py(idx(1)));;
        end
end
    
toc;

DATA = MakeFits(DATA,[]);
if DATA.showarea == 1
    DATA.lunatefit.fit= DATA.usefit.lunate;
    DATA.lunatefit.qfit= DATA.qfit.lunate;
end
set(DATA.top,'UserData',DATA);

DATA = RePlot(DATA);

function [ssq, p] = FitPlane(x, X,Y,Z)

id = find(~isnan(Z));
p = x(1) + x(2) .* X + x(3) .* Y;
ssq = sum((Z(id) - p(id)).^2);

function PlotGridImage(DATA, type)
[X,Y] =meshgrid(floor(min(DATA.px)):ceil(max(DATA.px)),floor(min(DATA.py)):ceil(max(DATA.py)));
rfxy(1,:) = DATA.map.rf(:,1)-DATA.map.rf(:,9);
rfxy(2,:) = DATA.map.rf(:,2)-DATA.map.rf(:,10);
rfxy(3,:) = DATA.map.rf(:,7);
rfxy(4,:) = DATA.map.rf(:,8);
for j = 1:size(X,1)
    for k = 1:size(X,2)
        id = find(abs(DATA.px - X(j,k)) < 0.7 & abs(DATA.py - Y(j,k)) < 0.7);
        px = DATA.px(id);
        py = DATA.py(id);
        if ~isempty(id)
            xid = find(DATA.map.pen(:,2) == X(j,k) & DATA.map.pen(:,3) == Y(j,k)) 
            age = min(DATA.map.age(xid)); %age is days back from now
            if type == 15 && ~isempty(xid)
                Z(j,k) = age;
            else
                Z(j,k) = mean(DATA.xc(id));
                ZY(j,k) = mean(DATA.yc(id));
            end
        else
            Z(j,k) = NaN;
            ZY(j,k) = NaN;
        end
    end
end
GetFigure('PenAge','parent',DATA.top);
[X,Y,Z] = fillpmesh(X,Y,Z);
hold off;
h = pcolor(X,Y,Z);
set(h,'buttondownfcn',@HitImage);
colorbar;




function DATA = PlotRFContour(DATA, type,varargin)
strargs = cell2cellstr(varargin);

showdiffs =0;
if sum(strcmp('diffs',strargs))
    showdiffs = 1;
end
[X,Y] =meshgrid(floor(min(DATA.px)):ceil(max(DATA.px)),floor(min(DATA.py)):ceil(max(DATA.py)));
rfxy(1,:) = DATA.map.rf(:,1)-DATA.map.rf(:,9);
rfxy(2,:) = DATA.map.rf(:,2)-DATA.map.rf(:,10);
rfxy(3,:) = DATA.map.rf(:,7);
rfxy(4,:) = DATA.map.rf(:,8);
for j = 1:size(X,1)
    for k = 1:size(X,2)
        id = find(abs(DATA.px - X(j,k)) < 0.7 & abs(DATA.py - Y(j,k)) < 0.7);
        px = DATA.px(id);
        py = DATA.py(id);
        if ~isempty(id)
            Z(j,k) = mean(DATA.xc(id));
            ZY(j,k) = mean(DATA.yc(id));
        else
            Z(j,k) = NaN;
            ZY(j,k) = NaN;
        end
        if type == 13 && ~isempty(id)
            idx = find(abs(rfxy(3,:)-px) < 0.5 & abs(rfxy(4,:)-py) < 0.5 & DATA.selected);
            if length(idx) > 1
                r = abs((rfxy(1,idx)-DATA.xc(j)) + i .* (rfxy(2,idx)-DATA.yc(j)));
                ZY(j,k) = std(r);
                if ZY(j,k) > 1
                    n = length(r);
                end
            else
                ZY(j,k) = NaN;
            end
        end
    end
end

if ismember(type,[5 14])
    tfit = [];
    yrange = [floor(min(min(Y))) ceil(max(max(Y)))];
    xrange = [floor(min(min(X))) ceil(max(max(X)))];
    [Xi, Yi] = meshgrid(xrange(1):0.1:xrange(2), yrange(1):0.1:yrange(2));
    if type == 5
        Zi = Interpf(X,Y, Z, Xi, Yi, 1, 0.5);
    elseif type == 14
        GetFigure('Fitted Map','parent',DATA.toplevel);
        hold off;
        if strcmp(DATA.fittype,'planar')
            [Z, ZY] = FitTopography(DATA.xyfit,X, Y);
            tfit = DATA.xyfit;
        elseif strcmp(DATA.fittype,'quadratic')
            [Z, ZY] = FitTopography(DATA.qfit,X, Y);
            tfit = DATA.qfit;
        elseif strcmp(DATA.fittype,'cubic')
            [Z, ZY] = FitTopography(DATA.cubicfit,X, Y);
            tfit = DATA.cubicfit;
        elseif strcmp(DATA.fittype,'rotate')
            [Z, ZY] = FitTopography(DATA.fit,X, Y);
            tfit = DATA.fit;
        else
            [Z, ZY] = FitTopography(DATA.cubicfit,X, Y);
            tfit = DATA.xyfit;
        end
    end
    [c, h] = contour(X,Y,Z,floor(min(min(Z))):ceil(max(max(Z))),'r');
    clabel(c,h,'color','r');
    hold on;
    [c, h] = contour(X,Y,ZY,floor(min(min(ZY))):ceil(max(max(ZY))),'b');
    clabel(c,h,'color','b');
    if isfield(tfit,'lunate')
        plot(tfit.lunate(1,2:end),tfit.lunate(2,2:end),'k-');
    end
    if showdiffs
        gid = find(DATA.selected > 0);
        rx = diff(DATA.map.rf(gid,[9 1]),[],2);
        ry = diff(DATA.map.rf(gid,[10 2]),[],2);
        id = find(~isnan(rx) & ~isnan(ry) & ~isnan(DATA.map.rf(gid,7)) & ~isnan(DATA.map.rf(gid,8)));
        if isfield(tfit,'excluded')
            id = setdiff(id,tfit.excluded);
        end
        rf = rx(id) + i * ry(id);
        d = abs(rf - mean(rf));
        crit = prctile(d,90).*3;
        gid = gid(id(d < crit));
        rx = rx(id(d<crit));
        ry = ry(id(d<crit));
        px = DATA.map.rf(gid,7);
        py = DATA.map.rf(gid,8);
        [Z, ZY] = FitTopography(tfit,rx,ry,'invert');
        d = (px + i*py) - (Z' + i*ZY');
        for j = 1:length(rx)
            plot([px(j) Z(j)], [py(j) ZY(j)],'-');
        end
        [a,b] = sort(abs(d));
    end
    if isfield(tfit,'lunate')
        DATA.lunatefit.fit = tfit.lunate;
    end
elseif type == 6
    GetFigure('XPcolor','parent',DATA.top);
    guess = [5 0.8 0.8];
    [ssq, fitz] = FitPlane(guess,X,Y,Z);
    options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
    fit = fminsearch(@FitPlane, guess, options, X,Y,Z);
    [ssq, fitz] = FitPlane(fit,X,Y,Z);
    DATA.fitmap.X = X;
    DATA.fitmap.Y = Y;
    DATA.fitmap.xposfit = fit;

    [X,Y,Z] = fillpmesh(X,Y,Z);
    hold off;
    h = pcolor(X,Y,Z);
    set(h,'buttondownfcn',@HitImage);
    colorbar;
elseif type == 7
    GetFigure('YPcolor','parent',DATA.top);
    DATA.fitmap.X = X;
    DATA.fitmap.Y = Y;
    guess = [2 0.3 -0.5];
    [ssq, fitz] = FitPlane(guess,X,Y,Z);
    options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
    fit = fminsearch(@FitPlane, guess, options, X,Y,Z);
    [ssq, fitz] = FitPlane(fit,X,Y,Z);
    DATA.fitmap.yposfit = fit;
    
    [X,Y,Z] = fillpmesh(X,Y,ZY);
    hold off;
    h = pcolor(X,Y,Z);
    set(h,'buttondownfcn',@HitImage);
    colorbar;
elseif type == 13
    GetFigure('RFscatter',DATA.top);
    hold off; 
    [X,Y,Z] = fillpmesh(X,Y,ZY)
    h = pcolor(X,Y,Z);
    set(h,'buttondownfcn',@HitImage);
    colorbar;
end

function HitImage(a,b)

DATA = GetDataFromFig(a);
c = get(gca,'CurrentPoint');
px = floor(c(1));
py = floor(c(1,2));
if isfield(DATA,'fit')
    [xy(1), xy(2)] = FitTopography(DATA.fit,px,py);
    str = sprintf(' RF from map fit %.1f,%.1f',xy);
else
    str = '';
end
fprintf('At %d,%d %s\n',px,py,str);
plotpen(DATA,px,py);


function PlotRFScatter(DATA)

map = DATA.map
pens = map.pen;
rfxy(1,:) = map.rf(:,1)-map.rf(:,9);
rfxy(2,:) = map.rf(:,2)-map.rf(:,10);
for j = 1:length(DATA.xc)
 px = DATA.px(j);
 py = DATA.py(j);
    idx = find(abs(pens(:,2) - px) < 0.6 & abs(pens(:,3) - py) < 0.6);

    if length(idx) >  1
    end
end

function PlotRFMap(DATA,gridlines)


MapDefs;
GetFigure(DATA.tags.fig);
hold off;
area = DATA.plot.area;
top = num2str(double(DATA.top));

if gridlines
%Draw blue lines for rows, and Red lines connecting columns
for x= unique(DATA.px)
    idx = find(DATA.px == x);
    plot(DATA.xc(idx),DATA.yc(idx),':');
    hold on;
end
for y= unique(DATA.py)
    idx = find(DATA.py == y);
    plot(DATA.xc(idx),DATA.yc(idx),'r:');
end
end

offset = range(DATA.xc)/60;

for j = 1:length(DATA.xc)
    if area == COMPAREV1V2
        plot([DATA.xc(j) DATA.axc(j)],[DATA.yc(j) DATA.ayc(j)],'-');
        hold on;
        plot(DATA.axc(j),DATA.ayc(j),'ro','buttondownfcn',['PlotMap(' top ',''point'',' num2str(j) ',1);']);
    end
    h = plot(DATA.xc(j),DATA.yc(j),'o','buttondownfcn',['PlotMap(' top ',''point'',' num2str(j) ',1);']);
    if DATA.penmarked(j)
        set(h, 'MarkerFaceColor','b');
    elseif DATA.ptype(j) == MULTICONTACT
        set(h, 'MarkerFaceColor','r');
    end
    
    if DATA.plot.labelpts
        if mod(DATA.px(j),1) > 0.1 | mod(DATA.py(j),1) > 0.1
            text(DATA.xc(j)+offset,DATA.yc(j),sprintf('%.1f,%.1f',DATA.px(j),DATA.py(j)),'fontsiz',DATA.plot.fontsiz);
        else
            text(DATA.xc(j)+offset,DATA.yc(j),sprintf('%.0f,%.0f',DATA.px(j),DATA.py(j)),'fontsiz',DATA.plot.fontsiz);
        end
    end
end

axis('equal');
if isfield(DATA.map,'drawlines')
    for j = 1:length(DATA.map.drawlines)
        plot(DATA.map.drawlines(j).x,DATA.map.drawlines(j).y,'-');
    end
end



function PlotGrid(DATA)

MapDefs;

top = num2str(double(DATA.top));


for j = 1:length(DATA.xc)
    if DATA.plot.labelpts
        text(DATA.px(j),DATA.py(j),sprintf('%.1f,%.1f',DATA.xc(j),DATA.yc(j)),'fontsiz',DATA.plot.fontsiz,'HorizontalAlignment','center','VerticalAlign','bottom');
    end
    h = plot(DATA.px(j),DATA.py(j),'bo','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(j) ',1);'],'MarkerFaceColor','b');
    if DATA.penmarked(j)
        set(h, 'MarkerFaceColor','m');
    elseif DATA.ptype(j) == MULTICONTACT
        set(h, 'MarkerFaceColor','r');
    end
    hold on;
end
%id is the penetrations selecte (e.g.. by visual area).
%DATA.pid is all penetration numbers
pid = find(ismember(DATA.id,DATA.pid));
[sorted, idx] = sort(DATA.id(pid));
idx = pid(idx);
colors = {[0 0 0], [0.3 0 0], [0.6 0.3 0] [0.6 0 0], [1 1 0 ], [1 0 0]} ;

k = 1;
if length(idx) > 5
    back = 5;
else
    back = length(idx)-1;
end
for j = idx(end-back:end)
    p = find(DATA.pid == DATA.id(j));
    hs(k) = plot(DATA.ppx(p),DATA.ppy(p),'o','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(DATA.ppid(p)) ',1);'],'color',colors{k},'MarkerFaceColor',colors{k});
    labels{k} = sprintf('%d',DATA.id(j));
    k = k+1;
end
legend(hs,labels);
for j =1:length(DATA.missed)
    k = find(DATA.map.pen(:,1) == DATA.missed(j));
    h = plot(DATA.map.pen(k,2),DATA.map.pen(k,3),'o','buttondownfcn',['PlotMap(' top ',''missedpoint'',' num2str(j) ',1);'],'color','k','MarkerFaceColor','none');
end

if isfield(DATA.map,'drawlines')
    for j = 1:length(DATA.map.drawlines)
        x = DATA.map.drawlines(j).x;
        y = DATA.map.drawlines(j).y;
        if length(x) < 100
            [x,y] = SampleLine(x,y);
        end
        [fx, fy] = FitTopography(DATA.xyfit,x,y,'invert');
        plot(fx,fy,'-');
    end
end
if isfield(DATA,'lunatefit')
    plot(DATA.lunatefit.fit(1,2:end),DATA.lunatefit.fit(2,2:end),'k-');
end

xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');
set(gca,'buttondownfcn',@HitImage);
setappdata(gcf,'ParentFigure',DATA.top);

function [sx,sy] = SampleLine(x,y, varargin)
npts = 10;
sx = [];
sy = [];
for j = 1:length(x)-1
    sx = [sx linspace(x(j),x(j+1),npts)];
    sy = [sy linspace(y(j),y(j+1),npts)];
end



function PlotGMDepths(DATA, varargin)

np = 0;
smooth = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'smooth',4)
        j = j+1;
        smooth = varargin{j};
    end
    j = j+1;
end
for j = 1:length(DATA.map.pens)
    pen = DATA.map.pens{j};
    if isfield(pen,'cmtype')
    cid = find(pen.cmtype == 6);
    mtcells = [];
    n = 0;
    for k = 1:length(cid)
        if strfind(pen.comments{cid(k)},'MT')
            n = n+1;
            mtcells(n) = pen.cmtime(cid(k));
        end
    end
    if n
        id = find(pen.cmtype == 5);
        gmid = find(pen.cmtime(id) <= mtcells(1));
        if length(gmid) && now - pen.datenum < DATA.plot.maxage
            gmid = gmid(end);
            if strncmpi(pen.comments{id(gmid)},'Grey',4)
            np = np+1;
            penx(np) = pen.pos(1);
            peny(np) = pen.pos(2);
            penz(np) = pen.depths(pen.cmtime(id(gmid))) - pen.enterdepth;
            pid(np) = j;
            else
                np = np+1;
                penx(np) = pen.pos(1);
                peny(np) = pen.pos(2);
                penz(np) = pen.depths(pen.cmtime(cid(1))) - pen.enterdepth;
                pid(np) = j;
            end
        end
    end
    end
end
if smooth > 0
[X,Y] = meshgrid(linspace(min(penx),max(penx)),linspace(min(peny),max(peny)));
Z = Interpf(penx,peny,penz,X,Y,1,smooth);
pcolor(X,Y,Z);
shading('interp');
colorbar;
else
    for j = 1:length(penz)
        plot3(penx(j),peny(j),penz(j),'o','buttondownfcn',['PlotMap(' num2str(DATA.top) ',''penpoint'',' num2str(pid(j)) ',1);']);
        hold on;
    end
end




function PlotGridArrows(DATA, varargin)

MapDefs;


plotmean = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'mean',4)
        plotmean = 1;
    end
    j = j+1;
end

colors = {[0 0 0], [0.3 0 0], [0.6 0.3 0] [0.6 0 0], [1 1 0 ], [1 0 0]} ;
top = num2str(DATA.top);
scale = 10;
id = find(DATA.selected);
for k = 1:length(id)
    j = id(k);
    p = DATA.map.pen(j);
    px = DATA.map.pen(j,2);
    py = DATA.map.pen(j,3);
    xs(:,k) = [px px+DATA.map.rf(j,1)./scale];
    ys(:,k) = [py py+DATA.map.rf(j,2)./scale];
    arrow([px px+DATA.map.rf(j,1)./scale],[py py+DATA.map.rf(j,2)./scale],20,0.1);
    hold on;
    pe = find(DATA.id == p);
    if length(pe) == 1
    ahs(k) = plot(px,py,'o', 'markerfacecolor','b','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(pe) ',1);']);
    end
    hs(k) = plot(px+DATA.map.rf(j,1)./scale,py+DATA.map.rf(j,2)./scale,'o','buttondownfcn',['PlotMap(' top ',''cellpoint'',' num2str(j) ',1);']);
    text(0.1+px+DATA.map.rf(j,1)./scale,py+DATA.map.rf(j,2)./scale,DATA.map.cellname{j});
end

if plotmean
    hold off;
    pens = xs(1,:) + i* ys(1,:);
    up = unique(pens);
    for j = 1:length(up)
        id = find(pens == up(j));
        xm = mean(xs(2,id));
        ym = mean(ys(2,id));
        px = real(up(j));
        py = imag(up(j));
        ahs(k) = plot(px,py,'o', 'markerfacecolor','b','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(pe) ',1);']);
        arrow([px xm],[py ym],20,0.1);
        hold on;
    end
end

minx = floor(min(xs(:)))-1;
maxx = ceil(max(xs(:)))+1;
miny = floor(min(ys(:)))-1;
maxy = ceil(max(ys(:)))+1;
set(gca,'xlim',[minx maxx]);
set(gca,'ylim',[miny maxy]);
[sorted, idx] = sort(DATA.pid);
colors = {[0 0 0], [0.3 0 0], [0.6 0.3 0] [0.6 0 0], [1 1 0 ], [1 0 0]} ;

if length(idx) > 5
    back = 5;
else
    back = length(idx)-1;
end
k = 1;
hs = [];
for j = idx(end-back:end)
    hs(k) = plot(DATA.ppx(j),DATA.ppy(j),'o','buttondownfcn',['PlotMap(' top ',''gridpoint'',' num2str(DATA.ppid(j)) ',1);'],'color',colors{k},'MarkerFaceColor',colors{k});
    labels{k} = sprintf('%d',DATA.pid(j));
    k = k+1;
end
legend(hs,labels);
if ~isempty(DATA.missed)
id = find(DATA.missed == MT);
for j =1:length(id)
    k = find(DATA.map.pen(:,1) == id(j));
    h = plot(DATA.map.pen(k,2),DATA.map.pen(k,3),'co','buttondownfcn',['PlotMap(' top ',''missedpoint'',' num2str(id(j)) ',1);'],'color','c','MarkerFaceColor','c');
end
id = find(DATA.missed == STS);
for j =1:length(id)
    k = find(DATA.map.pen(:,1) == id(j));
    h = plot(DATA.map.pen(k,2),DATA.map.pen(k,3),'co','buttondownfcn',['PlotMap(' top ',''missedpoint'',' num2str(id(j)) ',1);'],'color','g','MarkerFaceColor','g');
end
end


xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');


function PlotGridVar(DATA)

for j = 1:length(DATA.xc)
    if DATA.plot.labelpts
        text(DATA.px(j),DATA.py(j),sprintf('%.1f,%.1f',DATA.xc(j),DATA.xc(j)),'fontsiz',DATA.plot.fontsiz,'HorizontalAlignment','center','VerticalAlign','bottom');
    end
    h = plot3(DATA.px(j),DATA.py(j),DATA.rvar(j),'o','buttondownfcn',['PlotMap(DATA.top,''gridpoint'',' num2str(j) ',1);'],'MarkerFaceColor','b');
    if DATA.penmarked(j)
        set(h,'MarkerFaceColor','r');
    end
    hold on;
end
xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');

function PlotGridDepth(DATA)

hold off;
idx = find(DATA.selected > 0);
selected = idx;
if isempty(idx)
    fprintf('No data for Area %d\n',area);
    return;
end

for j = idx;
    if(DATA.map.depth(j) > -2000)
    if DATA.plot.labelpts
        text(DATA.px(j),DATA.py(j),sprintf('%.1f,%.1f',DATA.xc(j),DATA.xc(j)),'fontsiz',DATA.plot.fontsiz,'HorizontalAlignment','center','VerticalAlign','bottom');
    end
    h = plot3(DATA.map.rf(j,1),DATA.map.rf(j,2),DATA.map.depth(j),'o','buttondownfcn',['PlotMap(DATA.top,''depthpoint'',' num2str(j) ',1);'],'MarkerFaceColor','b');
    if DATA.marked(j)
        set(h,'MarkerFaceColor','r');
    end
    hold on;
    end
end
xr = get(gca,'Xlim');
yr = get(gca,'Ylim');
set(gca,'Xtick',[xr(1):xr(2)],'YTick',[yr(1):yr(2)],'Xgrid','on','Ygrid','on');

function ChoosePens(DATA)
  scrsz = get(0,'Screensize');
  wsc = scrsz(3)/1000;
  SPACE = 5;   VSPACE = 2;   BOXH = 15;
  
  selid = find(DATA.selected > 0)
  len = length(selid);
  if (len < 20)
      h = (len+1) * (15 + VSPACE);
  else
      h = (20+1) * (15 + VSPACE);
  end

  CloseTag('Blocks');
  cntrl_box = figure('Position', [200 scrsz(4)-(h+50) 500 h], 'Menubar', 'none',...
       'NumberTitle', 'off', 'Tag','Blocks','Name','BlockList');
  
   bp = [30 0 500 15];
   for k = 1:length(selid);
       j = selid(k);
      desc = sprintf('%d: %d,%d %s %.1f,%.1f',DATA.map.pen(j,1),DATA.map.pen(j,2),DATA.map.pen(j,3),...
          DATA.map.cellname{j},DATA.map.rf(j,2),DATA.map.rf(j,3));
      val = ~(DATA.excluded(j));
      uicontrol(gcf,'Style', 'CheckBox','String',desc,'Position', bp * wsc,...
      'Tag',sprintf('Block%d',k),'value',val,'UserData',j);
      bp(2) = bp(2) + bp(4) + VSPACE;
  end
  bp(3) = 25;
  bp(1) = 1;
  bp(2) = h-100;
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''checkblocks'')'],...
'String', 'Go', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);

  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''clearblocks'')'],...
'String', 'clr', 'Position', bp * wsc);
  bp(2) = bp(2)-bp(4);
  uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['PlotMap(' top ',''close'')'],...
'String', 'end', 'Position', bp * wsc);

  uicontrol(gcf, 'callback', ['PlotMap(' top ',''scroll'')'],'style','slider','min',1','max',k,'value',1,...
      'position',[480 5 20 h-10],'Tag','blockscroll');
  

function scrolllist(listpos)

BOXH = 15;
  VSPACE = 2;

blocks = GetBlockList('Block');
for j = 1:length(blocks)
    pos = get(blocks(j),'position');
    pos(2) = (j - (listpos)) * (BOXH + VSPACE);
    set(blocks(j),'position',pos);
end


function [blocks, isset, data] = GetBlockList(prefix, varargin)

j = 1;
it = 1;
while ~isempty(it)
it = findobj('Tag',sprintf('%s%d',prefix,j));
if ~isempty(it)
    blocks(j) = it;
    isset(j) = get(it,'value');
    data(j) = get(it,'UserData');
    j = j+1;
end
end

function ShowCell(map, id)
fprintf('%s %d %s %.1f %.1f %.1f\n',map.cellname{id},map.pen(id,1),map.datestr{id},map.rf(id,1),map.rf(id,2),map.depth(id));
  

function map = FixMap(map, fixes)


for j = 1:length(map.cellname)
    names{j} = GetName(map.cellname{j});
end

nl = 0;
for j = 1:length(fixes)
    name = GetName(fixes{j});
    id = find(strncmp(name,names,length(name)));
    for k = 1:length(id)
        if isfield(fixes{j},'pe')
            map.rf(id(k),6:8) = fixes{j}.pe(1:3);
            map.pen(id(k),1:3) = fixes{j}.pe(1:3);
        end
        if isfield(fixes{j},'area') 
            area = Name2Area(fixes{j}.area);
            if area > 0
                map.area(id(k)) = area;
            end
        end
    end
    if isfield(fixes{j},'line')
        nl = nl+1;
        map.drawlines(nl).x = fixes{j}.line(1,:);
        map.drawlines(nl).y = fixes{j}.line(2,:);
        map.drawlines(nl).type = 'line';
    end
    if isfield(fixes{j},'circle')
        nl = nl+1;
        map.drawlines(nl).xyr = fixes{j}.circle;
        map.drawlines(nl).type = 'circle';
        [map.drawlines(nl).x, map.drawlines(nl).y] = DrawEllipse(fixes{j}.circle,'noplot');
    end
end


function PrintPens(DATA,map)

[idx, order] = sort(map.datenum)
for j = 1:length(order)
    if DATA.selected(order(j))
    fprintf('%s P%d %.1f %.1f %s\n',map.datestr{order(j)},map.pen(order(j),1),map.pen(order(j),2),map.pen(order(j),3),map.cellname{order(j)});
    end
end

function SaveMap(map, monkey, varargin)

mapfile = sprintf('/b/data/%s/MapData.mat',monkey);
save(mapfile,'map');

function map = ReadAllPens(map, monkey, varargin)

pes = unique(map.pen(:,1));
reload = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'reload',4)
        reload = 1;
    end
    j = j+1;
end

penfile = sprintf('/b/data/%s/pens/pendata.mat',monkey);
penfile = CheckNameBug(penfile);
if exist(penfile,'file') && ~reload
    load(penfile);
    map.pens = pens;
    map.missed = missed;
    np = length(pens);
    for j = 1:length(pens)
        if isfield(pens{j},'num')
            oldpes(j) = pens{j}.num;
        end
    end
    id = find(oldpes < 0);
    oldpes(id) = 9000-oldpes(id); 
    id = find(pes < 0);
    pes(id) = 9000-pes(id); 
    newpes = setdiff(pes,oldpes);
else
    newpes = pes;
    np = 0;
end
for j = 1:length(newpes)
    pe = newpes(j);
    if pe >= 1
        name = sprintf('/b/data/%s/pens/pen%d.log',monkey,pe);
    else
        name = sprintf('/b/data/%s/pens/pen%d.log',monkey,9000-pe);
        pe = 9000-pe;
    end
    if isnan(pe)
        cprintf('red','Missing Pen number\n');
    else
    map.pens{pe} = ReadPen(name,'noplot');
    if ~isfield(map.pens{pe},'num')
        map.pens{pe}.num = pe;
    end
    if isfield(map.pens{pe},'missed')
        map.missed(pe) = map.pens{pe}.missed;
    else
        map.missed(pe) = 0;
    end
    end
end
if length(newpes)
    pens = map.pens;
    missed = map.missed;
    try
    save(penfile,'pens','missed');
    catch ME
        CheckExceptions(ME);
    end
end
pes = pes(pes > 0);
for j = 1:length(pes)
    pe = pes(j);
    if isfield(map.pens{pe},'boundary')
        B = map.pens{pe}.boundary;
        id = find(map.rf(:,6) == pe & map.depth' > B.depth);
        map.area(id) = Name2Area(B.below);
        id = find(map.rf(:,6) == pe & map.depth' < B.depth);
        map.area(id) = Name2Area(B.above);
    end
end