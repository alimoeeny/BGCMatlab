function PlotDuplicates(DATA, varargin)

A = getappdata(DATA.toplevel,'AllDuplicates');

PC.PlotCellList(DATA,'nomarks','calltag','duplicates');

Im = [];
colors = mycolors;
for j = 1:length(A)
    if isfield(A{j},'groups')
    G = A{j}.groups;
    nd = 0;
    for k = 1:length(G)
        if length(G{k}) > 1
            nd = nd+1;
            cells = A{j}.cellid(G{k});
            p = A{j}.probes(G{k});
            Im(j,floor(p)) = k;
            for b = 1:length(G{k})
                PC.DrawBox(j,p(b),'duplicates','color',colors{nd},'linewidth',2);
            end
        end
    end
    end
end
GetFigure('Duplcates','parentfigure',DATA.toplevel);
h = imagesc(Im);
set(h,'buttondownfcn',{@PC.HitImage, 'duplicates'});



function HitImage(src, b, varargin)

DATA = GetDataFromFig(src);
A = getappdata(DATA.toplevel,'AllDuplicates');
ax = get(src,'Parent');
xy = get(ax,'currentpoint');
l = get(ax,'Children');
tag = get(get(ax,'Parent'),'Tag');
ex = round(xy(1,2));
p = round(xy(1,1));
if xy(1,1) > p
    setcl = 2;
else
    setcl = 1;
end
[Clusters, DATA] = PC.CheckClusterLoaded(DATA, ex);
