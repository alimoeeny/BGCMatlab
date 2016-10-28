function DATA = SaveLayout(DATA, varargin)
%DATA = SaveLayout(DATA, varargin)
%Saves Layout of windows for application
%with tags labelled in DATA.tag
%uses DATA.layoutfile, or file named with
%SaveLayout(DATA, 'layout', file
%       if file is not full pathname, then looks in
%        DATA.prefsdir (if exists)
%SaveLayout(DATA, 'nochoose')  Does not bring up a gui for selecting
%finename
% see also ApplyLayout

choosefile = 1;
layoutfile = [];
if isfield(DATA,'layoutfile')
    layoutfile = DATA.layoutfile;
elseif isfield(DATA,'defaultlayout')
    layoutfile = DATA.defaultlayout;    
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'layout',5) || strncmpi(varargin{j},'name',4)
        j = j+1;
        layoutfile = varargin{j};
    elseif strncmpi(varargin{j},'folder',5) %build name with ui and hostname
        j = 1+1;
        savedir = varargin{j};
        layoutfile = [savedir '/' gethostname '.' GetUserName];
    elseif strncmpi(varargin{j},'nochoose',5)
        j = j+1;
        choosefile = 0;
    elseif j == 1 && ischar(varargin{j})
        layoutfile = varargin{j};
    end
    j = j+1;
end


savefields = {'Figpos'};
if isfield(DATA,'gui')
    GUI = DATA.gui;
    savefields = {savefields{:} 'GUI'};
end
proceed = 1;
Figpos = GetFigPos(DATA.toplevel); %sets the app data, appending if necessary
Figpos = getappdata(DATA.toplevel,'Figpos'); 
if isfield(DATA,'tag')
f = fields(DATA.tag);
   for j = 1:length(f);
       it = findobj('type','figure','Tag', DATA.tag.(f{j}));
       if length(it) == 1
           Figpos.(DATA.tag.(f{j})) = get(it,'Position');
       end
   end
end
   if choosefile
       [outname, pathname] = uiputfile(layoutfile);
       if outname
           DATA.layoutfile = [pathname outname];
       else
           proceed = 0;
       end
   else
       DATA.layoutfile = layoutfile;
   end
   if proceed
       try
       save(CheckNameBug(DATA.layoutfile),savefields{:});
       fprintf('Layout saved to %s\n',DATA.layoutfile);
       catch
           cprintf('red','Error Saving  %s\n',DATA.layoutfile);
       end
   end
   setappdata(DATA.toplevel,'Figpos',Figpos);
   