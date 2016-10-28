function DATA = SaveCallback(DATA, a, varargin)
%DATA = SaveCallback(DATA, a, varargin)
%see also PopupWindow;
if ishandle(a)
    if isfield(DATA,'interactive') && DATA.interactive == 0
        DATA.interactive = 1;
    end
    if isfield(DATA,'guistate')
        DATA.guistate.lastguihandle = a;
    else
        DATA.gui.lastguihandle = a;
    end
    SetData(DATA);
end