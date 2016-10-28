function [ex,p,setcl] = Mouse2Expt(src, data)
%PC.Mouse2Expt convert mouse click to eid, pid.

    ax = gca;
xy = get(ax,'currentpoint');
srctype = get(src,'type');
srcdata = get(src,'UserData');
l = get(ax,'Children');
tag = get(get(ax,'Parent'),'Tag');
ex = round(xy(1,2));
p = round(xy(1,1));
if xy(1,1) > p
    setcl = 2;
else
    setcl = 1;
end