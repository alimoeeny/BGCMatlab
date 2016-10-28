function xc = BuildShapeMap(DATA, C, Clusters, celln, varargin)
%PC.BuildShapeMap(DATA, CellList, Clusters, cell) calc Xcorrs for
%all probes and expts for a given cell, then plot map


[e,p,cl] = celllist.find(C,celln);
[a,b] = Counts(p);
[c,d] = max(a);
bestp = b(d);
id = find(p == (b(d))); 
A = PC.GetClusterInfo(Clusters,[e(id(1)) p(id(1)) cl(id(1))]);
ms = A.MeanSpike.ms .* A.ncut;
ns = A.ncut;
for j = 2:length(id)
    A = PC.GetClusterInfo(Clusters,[e(id(1)) p(id(1)) cl(id(1))]);
    ms = ms + A.MeanSpike.ms .* A.ncut;
    ns = ns + A.ncut;
end
ms = ms./ns;
im = [];
A.MeanSpike.ms = ms;

for e = 1:size(C,1)
    for p = 1:size(C,2)
        for c = 1:size(C,3) %all ecps
            B = PC.GetClusterInfo(Clusters,e,p,c,'allexpt');
            [pxc(c), details] = PC.ShapeCorr(A,B,'shifts',5,'drifts',5,'proberange',5);
            ppx(c) = details.probeshift;
        end
        [xc(e,p), clid(e,p)] = max(pxc);
        pshift(e,p) = ppx(clid(e,p));
    end
end

GetFigure('CompareShape','parent',DATA.toplevel);
h = imagesc(xc);
colorbar;
caxis([0.5 1]);
setappdata(gcf,'CellMeanSpike',A);
setappdata(gcf,'BestClusterid',clid);
set(h,'buttondownfcn',{@HitImage});

function HitImage(src,b)

DATA = GetDataFromFig(src);
xy = get(gca,'currentpoint');
ex = round(xy(1,2));
p = round(xy(1,1));

clid = getappdata(gcf,'BestClusterid');
fprintf('Comparing E%dP%dcl%d\n',ex,p,clid(ex,p));
A = getappdata(gcf,'CellMeanSpike');
Clusters = getappdata(DATA.toplevel,'AutoClusters');
B = PC.GetClusterInfo(Clusters, [ex p clid(ex,p)]);
PC.CompareSpikeShape({A B},'tag','CellMeanShape');

