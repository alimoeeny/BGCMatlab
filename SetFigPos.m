function SetFigPos(Figpos, tag, varargin)
%SetFigPos(X, tag)
%Sets figure size and location for figure whose tag is tag
%Data used to determine figure pos depends on X
%If X is a figure, getappdata(X,'Figpos');
%If X is strucutre with a field name matching tag, then 
%thie field of the structure is used
%if X is a typical figure data structure withe the field 'toplevel,
%then getappdata(X.toplevel,'Figpos');
%SetFigPos(0) restores all figure locations last stored with GetFigPos(0);
% See also GetFigPos

strargs = cell2cellstr(varargin);
toplevel = 0;
apptag = [];
if sum(strcmp('force',varargin))
    forcepos = 1;
else
    forcepos = 0;
end
newtag = '';
if isfigure(Figpos)
    if isappdata(Figpos,'ParentFigure')
        toplevel = getappdata(Figpos,'ParentFigure');
        if nargin >1 && strcmp(tag,'force')
            forcepos = 1;
        end
        newtag = get(Figpos,'tag');
    else
        toplevel = Figpos;
    end
    Figpos = getappdata(toplevel,'Figpos');
elseif isfield(Figpos,'toplevel')
    toplevel = Figpos.toplevel;
    Figpos = getappdata(Figpos.toplevel,'Figpos');
elseif isnumeric(Figpos) && Figpos == 0
    toplevel = 0;
    Figpos = getappdata(0,'Figpos');
end
if nargin == 1
    f = fields(Figpos)
   for j = 1:length(f)
       F = findobj(allchild(0),'flat','type','figure','tag',f{j});
       if ~isempty(F)
           set(F,'position',Figpos.(f{j})(1:4));
       elseif strncmp(f{j},'Figure',6)
           F = sscanf(f{j},'Figure%d');
           set(F,'position',Figpos.(f{j})(1:4));
       end
   end
   return;
end
    it = findobj('type','figure','Tag',tag);
    if isempty(it)       
        it = findobj('type','figure','Tag',newtag);
        if isempty(it)
            return;
        end
    end
    
if isappdata(toplevel,'ApplicationTag')
    apptag = getappdata(toplevel,'ApplicationTag');
end
if isfield(apptag,'extra')
    tag = strrep(tag,apptag.extra,'');
end
f = genvarname(tag);


if isfield(Figpos,f)
    go = 1;
    if length(Figpos.(f)) > 4 && Figpos.(f)(5) == 1 %already set once - user mmay have moved
        go = forcepos;
    end
    if length(it) == 1 && go
        set(it,'position',Figpos.(f)(1:4));
        Figpos.(f)(5) = 1;
    end
else
   Figpos.(f)  = get(it(1),'position');
end
if isfigure(toplevel)
    setappdata(toplevel,'Figpos',Figpos);
end
