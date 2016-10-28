function PopupWindow(toplevel, tag, label, callback, varargin)

DATA = get(toplevel,'UserData');
if isfield(DATA,'tag') && isfield(DATA.tag,'popup')
    wtag = DATA.tag.popup;
else
    wtag = 'popup';
end

if nargin == 1
    tag = '';
    label = '';
    callback = [];
end

truetag = tag;
truelabel = label;
if isempty(tag)
    tag = label;
end
tag = strrep(tag,'&','');
label = strrep(label,'&','');
X = getappdata(toplevel,'PopupWindowData');
if ~iscell(X) 
    X = {};
    if nargin ==1 % dont popup if empty
    return;
    end
end
[F, isnew] = GetFigure(wtag,'parent',toplevel,'trackpos','setpos');
if isnew
    X = RemoveDuplicates(X);
    set(F,'menubar','none');
    rm = uimenu(F,'label','Remove','tag','RemoveMenu');
    sm = uimenu(F,'label','Add','tag','AddMenu');
    uimenu(F,'label','Dock','tag','DockMenu','callback',@DockButtons);
    uimenu(sm,'label','Last Operation',...
        'callback',@AddCallback,...
        'tag','AddCallbackButton');
    if size(X,2) > 2
        a = CellToStruct(X(:,3));
        nr = 1+sum([a.on]);
        id = find([a.on]);
        if sum(strcmp(truelabel,X(id,2))) == 0
            nr = nr-1;
        end
    else
        nr = 1;
    end
    nb =0;
    for j = 1:size(X,1)
        X{j,1} = strrep(X{j,1},'&','');
        if X{j,3}.on
            it =FindChild(DATA.toplevel,'tag',X{j,1},'label',X{j,2});
            if ~isempty(it)
                X{j,3}.mid = it(1);
                nb = nb+1;
            else
                X{j,3}.mid = 0;
            end
        end
    end
    nr = nb;
    nb = 0;
    for j = 1:size(X,1)
        if X{j,3}.on
%make empty tags = label for new buttons            
            if isempty(X{j,1})
                X{j,1} = X{j,2};
            end
            if double(X{j,3}.mid) > 0
                nb = nb+1;
                X{j,3}.buttonid = uicontrol(F,'style','pushbutton','string',X{j,2},...
                    'callback',get(X{j,3}.mid,'Callback'),...
                    'tag',X{j,1},...
                    'units','normalized','position',[0.02 (nb-1)/nr 0.96 0.96/nr]);
                uimenu(rm,'label',X{j,2},'callback',{@RemoveButton, X{j,1}},'tag',X{j,1});
            end
            
        end
    end
end

it = findobj(F,'type','uicontrol','style','pushbutton','tag',tag);
if isempty(it) && ~isempty(tag)
    it = findobj(F,'type','uicontrol','style','pushbutton');
    nr = length(it)+1;
    for j = 1:length(it)
        set(it(j), 'units','normalized','position',[0.02 (j-1)/nr, 0.96, 0.96/nr]);
    end
    
    j = nr;
    if ~isempty(callback)
        state.on = 1;
    X{end+1,1} = truetag;
    X{end,2} = truelabel;
    state.mid = findobj(DATA.toplevel,'tag',tag');
    if length(state.mid) > 1
        state.mid = state.mid(1);
    end
    if strcmp(get(state.mid,'type'),'uimenu')
        state.menutype = 1;
    else
        state.menutype = 0; 
    end
    state.buttonid = uicontrol(F,'style','pushbutton','string',tag,...
        'callback',callback,...
        'tag',tag,...
        'units','normalized','position',[0.02 (j-1)/nr 0.96 0.96/nr],...
        'userdata',state);
    X{end,3} = state;
    end

    sm  = findobj(F,'type','uimenu','tag','RemoveMenu');
    uimenu(sm,'label',label,'tag',tag,'callback',{@RemoveButton, tag});
    setappdata(toplevel,'PopupWindowData',X);
end

function RemoveButton(a,b,tag)

F = GetFigure(a);
it = findobj(F,'type','uicontrol','style','pushbutton','tag',tag);
delete(it);
sm  = findobj(F,'type','uimenu','tag','RemoveMenu');
it = findobj(sm,'type','uimenu','tag',tag);
sm  = findobj(F,'type','uimenu','tag','AddMenu');
ait = findobj(sm,'type','uimenu','tag',tag);
if isempty(ait)
    args = CopyUiProperties(it(1),{'tag' 'label'});
    hm = uimenu(sm,args{:});
    tag = get(it,'tag');
    set(hm, 'callback', {@AddButton, tag},'userdata',get(it,'callback'));
end
delete(it);
parent = getappdata(F,'ParentFigure');
X = getappdata(parent,'PopupWindowData');
state = CellToStruct(X(:,3));
%[state.buttonid] can throw error if some are double, some handles
for j = 1:length(state)
  if state(j).buttonid == it
      X{j,3}.on = 0;
  end
end
setappdata(parent,'PopupWindowData',X);

function AddButton(a,b,tag)

F = GetFigure(a);
toplevel = getappdata(F,'ParentFigure');
sm = get(a,'parent');
%original callback stored in userdataeta
AllV.PopupWindow(toplevel, get(a,'tag'),get(a,'label'),get(a,'userdata'));


function AddCallback(a,b)

F = GetFigure(a);
DATA = GetDataFromFig(F);
a = DATA.guistate.lastguihandle;
if ishandle(a)
    AllV.PopupWindow(DATA.toplevel, get(a,'tag'),get(a,'label'),get(a,'callback'));
end

function DockButtons(a,b)
F = GetFigure(a);
DATA = GetDataFromFig(F);
X = getappdata(DATA.toplevel,'DockButtonData');
%remove any currently docked first
for j = 1:size(X,1)
    if isfield(X{j,3},'buttonid') && myhandle(X{j,3}.buttonid)
        delete(X{j,3}.buttonid);
    end
end
X = getappdata(DATA.toplevel,'PopupWindowData');
nc = size(X,1);
for j = 1:size(X,1)
    if X{j,3}.on
        it =FindChild(DATA.toplevel,'tag',X{j,1},'label',X{j,2});
        if ~isempty(it)
            state.mid = it(1);
            if strcmp(get(state.mid,'type'),'uimenu')
                state.menutype = 1;
            else
                state.menutype = 0;
            end
            X{j,3}.buttonid = uicontrol(DATA.toplevel,'style','pushbutton','string',X{j,2},...
                'callback',get(it(1),'Callback'),...
                'tag',X{j,1},...
                'units','normalized','position',[(j-1)/nc 0.0 0.96/nc 0.05],...
                'userdata',state);
        else
            X{j,3}.buttonid = 0; %no new button made
        end        
    end
end
setappdata(DATA.toplevel,'DockButtonData',X);

function X = RemoveDuplicates(X)

good = ones(size(X,1),1);
for j = 1:size(X,1)
    for k = 1:j-1
        if strcmp(X{j,1}, X{k,1}) && strcmp(X{j,2}, X{k,2})
            good(k) = 0;
        end
    end
end
good = find(good);
X = X(good,:);