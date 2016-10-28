function X = PopupWindow(a,b, op, varargin)

DATA = GetDataFromFig(a);
X = [];
ncol = 1;

popupmenu = { ...
    'cmpautomahal' 'Auto: compare mahal' {@PC.AutoCompare, 'mahal'}; 
    'Combine' 'Combine' {@PC.PlotMenu, 'expt', 'combine'}; 
    'Combinetest' 'Test Combine' {@PC.OptionMenu, 'combinexpt'}; 
    'RateSeq' 'Rate Sequnence' {@PC.PlotMenu, 'cells', 'rateseqone'}; 
    'autocompare','Compare All Auto' {@PC.OptionMenu, 'compareauto'};
    'cmpautosu' 'Auto: compare sus' {@PC.AutoCompare, 'sus'};
    'rateseqallcells' 'Rate Sequence All Cells',{@PC.PlotMenu, 'cells', 'rateseqall'};
    'checkrateseqs' 'Check &Rate sequences',{@PC.OptionMenu 'checkratesequence'};
    'excludecl1' 'Exclude Selected Trials',{@PC.OptionMenu, 'excludecl1'};
     'plotallexptprobe' 'Plot Expts for 1 Probe' {@PC.OptionMenu, 'setallprobeplot'};
     'plotallexptcell' 'Plot Expts for 1 Cell' {@PC.OptionMenu, 'setallcellplot'};
    'allexptspks' 'All Spks This Expt' {@PC.PlotMenu, 'probes', 'allspks'};
    'allprobespks' 'All Spks This Probe' {@PC.PlotMenu, 'probes', 'allprobespks'};
    'allexptxy' 'All &XY This Expt' {@PC.PlotMenu, 'probes', 'AllXY'};
    'allprobexy' 'All X&Y This Probe' {@PC.PlotMenu, 'probes', 'AllprobeXY'};
    'allprobemean' 'Mean Spikes' {@PC.PlotMenu, 'probes', 'AllMean'};
    'allcellxy' 'All Cells XY' {@PC.PlotMenu, 'cells', 'AllXY'};
    'allcellmean' 'All cells Mean' {@PC.PlotMenu, 'cells', 'AllMean'};
    'allcellmeanim' 'All Cells Mean Image' {@PC.PlotMenu, 'cells', 'AllMeanIm'};
    'allcellspks' 'AllSpks Cell' {@PC.PlotMenu, 'cells', 'allspks'};
    'allexptcellsspks' 'All Cells Spks Expt',{@PC.PlotMenu, 'probes', 'cellspks'};
    };

if strcmp(op,'addbutton')
    if isfield(DATA.gui,'xbutton')
    end
    x = [];
    if isfield(DATA.gui,'lastguihandle')
        x = DATA.gui.lastguihandle;
    end
    if ishandle(x)
    popupmenu = getappdata(DATA.toplevel,'popupmenu');
    blabel = regexprep(get(x,'label'),'\s','');
    if sum(strcmp(blabel,popupmenu(:,1))) == 0 % new button
        popupmenu{end+1,1} = blabel;
        popupmenu{end,2} = get(x,'label');
        popupmenu{end,3} = get(x,'callback');
        setappdata(DATA.toplevel,'popupmenu',popupmenu);
    end
    DATA.popup.(blabel) = 1;
    F = GetFigure(a);
    MakePopupButtons(DATA, F);
    SetData(DATA);
    end
elseif strcmp(op,'buildstruct');
    for j = 1:size(popupmenu,1)
        X.(popupmenu{j,1}) = 0;
    end
elseif strcmp(op,'rebuild')
    op = PC.PopupWindow(DATA,[],'buildstruct');
    DATA.popup = CopyFields(DATA.popup, op,'-ifnew');
    SetData(DATA);
elseif strcmp(op,'popup')
    [F, isnew] = GetFigure(DATA.tag.popup,'parent',DATA.toplevel,'trackpos');
    set(F,'menubar','none');
    if isnew
        setappdata(DATA.toplevel,'popupmenu',popupmenu);
        set(F,'menubar','none');
        hm = uimenu(F,'Label','Select');
        uimenu(hm,'Label','Check box list','Callback',{@PC.PopupWindow, 'select'});
        uimenu(hm,'Label','Add Button For Last Operation','Callback',{@PC.PopupWindow, 'addbutton'});
        uimenu(F,'Label','XX','Callback',{@PC.PopupWindow, 'rebuild'});
        PC.PopupWindow(DATA, [], 'select');
    end
    MakePopupButtons(DATA,F);
elseif strcmp(op,'setbutton')
    onoff = get(a,'value');
    DATA.popup.(popupmenu{varargin{1},1}) = onoff;
    PC.PopupWindow(DATA, [], 'popup');
    SetData(DATA);
elseif strcmp(op,'select')
    [F, isnew] = GetFigure([DATA.tag.popup 'select'],'parent',DATA.toplevel,'trackpos');
    if isnew
        set(F,'menubar','none');
    end
    MakePopCheckBox(DATA, F);
end


function MakePopupButtons(DATA,F)
    f = fields(DATA.popup); 
    nr = 0;
    popupmenu = getappdata(DATA.toplevel,'popupmenu');
    for j = 1:length(f)
        nr = nr+DATA.popup.(f{j});
    end
    bp(1) = 0.01;
    bp(2) = 0.01;
    bp(3) = 0.95;
    row = 0;
    for j = 1:length(f)
        id = find(strcmp(f{j},popupmenu(:,1)));
        if isempty(id)
            acknowledge(sprintf('Dont know how make button for\n%s\nAdd it Manually to Shortcut Popup',f{j}),'PlotClusters' )
        else
            if DATA.popup.(f{j})
                row = row+1;
                bp(4) = 1./nr;
                bp(2) =0.01 + (row-1)/nr;
                uicontrol(F,'style','pushbutton','String',popupmenu{id,2},'Callback',popupmenu{id,3},...
                    'units','normalized','position',bp);
            else
                it = findobj(F,'style','pushbutton', 'String',popupmenu{id,2});
                delete(it);
            end
        end
    end


function MakePopCheckBox(DATA, F)
    popupmenu = getappdata(DATA.toplevel,'popupmenu');
    nr = size(popupmenu,1);
    bp(1) = 0.01;
    bp(2) = 0.01;
    bp(3) = 0.95;
    bp(4) = 1./nr;
    for j = 1:size(popupmenu,1)
        bp(2) =0.01 + (j-1)/nr;
        uicontrol(F,'style','checkbox','String',popupmenu{j,2},'Callback',{@PC.PopupWindow, 'setbutton', j},...
            'value',DATA.popup.(popupmenu{j,1}),...
            'units','normalized','position',bp);
    end
