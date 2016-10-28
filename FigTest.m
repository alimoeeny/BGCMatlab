function F = FigTest(varargin)
%FigTest illustrates new behaviour of UImenus
%if a menu in the menubar has no submenus, it does not inactivate
%after a button press, but stays highlighted.  Just moving the mouse
%over another uimenu then calls the callback of THAT menu

F = figure;
h = uimenu(F,'label','Test','callback',@TestMark,'hittest','off');
h = uimenu(F,'label','Test2','callback',@TestMarkB,'hittest','off');
plot(cos(0:0.1:pi));
text(0,0.5,'Click On Test Menu, then move mouse over Test2 without clicking');

function TestMark(src, event)
fprintf('Hit1\n');

function TestMarkB(src, event)
fprintf('Hit2\n');