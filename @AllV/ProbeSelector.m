function ProbeSelector(DATA, type)if nargin == 1    type = 'set';end[F, isnew] = GetFigure(DATA.tag.probeselect,'parent',DATA.toplevel);   Clusters = AllV.mygetappdata(DATA,'Clusters');   marked = CellToMat(Clusters,'marked');   if length(marked) < DATA.allnprobes       marked(DATA.allnprobes) = 0;   end   if isfield(DATA.ArrayConfig,'badprobes')       marked(find(DATA.ArrayConfig.badprobes)) = 100;   end      if isfield(DATA,'CellList')   C = DATA.CellList;   C(C < 0) = 0;   C(C > 0) = 1;   C = sum(C,3);   nex = size(C,1);   C = sum(C,1)./nex;   else       C = zeros(1,DATA.allnprobes);   end   if sum(C>0.1) > 2  %more than 2 cells set        fc = [1 1 1];       bc = [0 0 0];   else       fc = [0 0 0];       bc = [1 1 1];   end   if isnew       nr = 10;       nc = 10;       w = 1./nc;       h = 1./nr;       X = uibuttongroup(F,'SelectionChangeFcn',{@AllV.SelectProbe, type},'tag','ProbeSelectorGroup');       set(X,'backgroundcolor',bc);       set(X,'foregroundcolor',fc);       set(X,'UserData',DATA.toplevel);       cmenu = uicontextmenu;       uimenu(cmenu,'label','Cell','foregroundcolor','y');       uimenu(cmenu,'label','Cell in other expts','foregroundcolor','c');       uimenu(cmenu,'label','Marked Good','foregroundcolor','r');       uimenu(cmenu,'label','Marked Good MU','foregroundcolor','b');       uimenu(cmenu,'label','Marked Bad/duplicate','foregroundcolor','g');       for j = 1:DATA.allnprobes           x = mod((j-1),nr)./nr;           y = 1 -1./nc - floor((j-1)./nc)/nc;           it = uicontrol(X,'style','radiobutton','string',num2str(j),'UserData', j,...               'units','norm','position',[x, y, w, h]);           [a, b] = AllV.isacell(DATA, DATA.exptno, j);           if a               set(it,'foregroundcolor','y');           elseif C(j) > 0 %Cells one other expts, this probe               set(it,'foregroundcolor','c');           elseif marked(j) == 2 %good               set(it,'foregroundcolor','r');           elseif marked(j) == 4 %good MU               set(it,'foregroundcolor','b');           elseif marked(j) == 3 || marked(j) == 5 || marked(j) == 100 %BAD               set(it,'foregroundcolor','g');           else               set(it,'foregroundcolor',fc);           end           set(it,'uicontextmenu',cmenu,'backgroundcolor',bc);       end       pos = get(it,'position');       bp(1) = pos(1)+pos(3);       bp(2) = pos(2);       bp(4) = pos(4);       bp(3) = 1-bp(1);       it = uicontrol(X,'style','text','string','Rclick for colorscheme','units','norm','position',bp,...           'foregroundcolor',fc,'backgroundcolor',bc);              set(it,'uicontextmenu',cmenu);       set(X,'uicontextmenu',cmenu);       set(gcf,'uicontextmenu',cmenu);   else       for j = 1:DATA.allnprobes           h = findobj(F,'style','radiobutton','UserData',j);           [a, b] = AllV.isacell(DATA, DATA.exptno, j);           if b               set(h,'foregroundcolor','y');           elseif C(j) > 0 %Cells one other expts, this probe               set(h,'foregroundcolor','c');           elseif ~isempty(h) && marked(j) == 2               set(h,'foregroundcolor','r');           elseif marked(j) == 4 %good MU               set(h,'foregroundcolor','b');           elseif marked(j) == 3 %BAD               set(h,'foregroundcolor','g');           else               set(h,'foregroundcolor',fc);           end       end   end        