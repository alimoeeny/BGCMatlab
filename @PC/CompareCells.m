function CompareCells(DATA, Clusters,C, cella, cellb)
%PC.CompareCells(DATA, Clusters,C, cella, cellb) check if cells are the
%smae
if isfigure(DATA)
    DATA = get(DATA,'UserData');
end
if isempty(Clusters)
    Clusters = getappdata(DATA.toplevel,'Clusters');
end

if size(cella,2) ==3 && size(cella,1) > 1 %list of points
    b = cellb;
    ae = cella(1:b,1);
    ap = cella(1:b,2);
    ac = cella(1:b,3);
    be = cella(2:end,1);
    bp = cella(2:end,2);
    bc = cella(2:end,3);    
else
    id = find(C == cella);
    [ae, ap, ac] = ind2sub(size(C),id);
    [ae,id] = sort(ae);
    ap = ap(id);
    ac = ac(id);
    id = find(C == cellb);
    [be, bp, bc] = ind2sub(size(C),id);
    [be,id] = sort(be);
    bp = bp(id);
    bc = bc(id);
    if isempty(ae)
        fprintf('Cell %d is empty\n',cella);
        return;
    end
    if isempty(be)
        fprintf('Cell %d is empty\n',cellb);
        return;
    end
    X.cells = [cella cellb];
end

    for j = 1:length(ae)
        for k = 1:length(be)
            A = PC.GetClusterInfo(Clusters,ae(j),ap(j),ac(j),'allexpt');
            B = PC.GetClusterInfo(Clusters,be(k),bp(k),bc(k),'allexpt');
            [xc, details{1}] = PC.ShapeCorr(A,B,'drift',1,'delays',4);
            M(j,k) = xc;
            exdiff(j,k) = abs(be(k)-ae(j));
%            M(k,j) = xc;
        end
    end
    [F, isnew] = GetFigure('CompareShape','parent',DATA.toplevel);
    if isnew
        delete(findobj(allchild(5),'flat','type','uimenu'));
        hm = uimenu(F,'label','Find');
        uimenu(hm,'label','Possible Pairs','callback',@FindPairs);
        hm = uimenu(F,'label',sprintf('Merge Cells %d %d',cellb,cella),'tag','cellmerge', 'callback',   {@CombineCells, cella, cellb});
    end
    it = findobj(allchild(F),'flat','tag','cellmerge');
    if ~isempty(it)
        set(it,'label',sprintf('Merge Cells%d %d',cellb,cella),'callback',   {@CombineCells, cella, cellb});
    end
    hold off;
    if length(ae) == 1 && length(be) ==1
        subplot(1,1,1);
        CompareSpikeShape(A,B);
    else
    subplot(1,2,1)        
    h = imagesc(M);
    set(h,'buttondownfcn',@HitImage);
    xlabel(sprintf('Expts for cell %d',cellb));
    ylabel(sprintf('Expts for cell %d',cella));
    xtick = 1:2:size(M,2);
    ytick = 1:2:size(M,1);
    set(gca,'xtick',xtick,'ytick',ytick);
    set(gca,'xticklabel',cellstr(num2str(DATA.CellDetails.exptids(be(xtick))')))
    set(gca,'yticklabel',cellstr(num2str(DATA.CellDetails.exptids(ae(ytick))')))
    colorbar;
    subplot(1,2,2);
    plot(exdiff(:),M(:),'o');
    X.alist = [ae ap ac];
    X.blist = [be bp bc];
    X.parent = DATA.toplevel;
    setappdata(gcf,'CellLists',X);
    end
    
    function CompareSpikeShape(A,B, varargin)
        c = mycolors;
        hold off;
        Sa = A.MeanSpike.ms;
        Sb = B.MeanSpike.ms;
        for j = 1:size(Sa,1)
            plot(Sa(j,:),'-','color',c{j});
            hold on;
            plot(Sb(j,:),'--','color',c{j});
        end
        [xca, details{1}] = PC.ShapeCorr(A,B,'drift',1);
        [xcb, details{1}] = PC.ShapeCorr(A,B,'delay',5,'drift',1);
        title(sprintf('Corr %.2f,%.2f',xcb,details{1}.truexc));
        if isfield(A,'cellnumber')
            cella = A.cellnumber;
        else
            cella = 0;
        end
        if isfield(B,'cellnumber')
            cellb = B.cellnumber;
        else
            cellb = 0;
        end
        labels{1} = sprintf('Cell%d E%dP%dcl%d',cella,A.exptno,A.probe(1),A.cluster);
        labels{2} = sprintf('Cell%d E%dP%dcl%d',cellb,B.exptno,B.probe(1),B.cluster);
        legend(labels);
        
function HitImage(a,b)
    F  = gcf;
    pos = get(gca,'currentpoint');
    xy = round([pos(1,1) pos(1,2)]);
    xy(xy<1) = 1;
    X = getappdata(F,'CellLists');
    ae = X.blist(xy(1),:);
    be = X.alist(xy(2),:);
    Clusters = getappdata(X.parent,'Clusters');
    A = PC.GetClusterInfo(Clusters,ae);
    A.cellnumber = X.cells(1);
    B = PC.GetClusterInfo(Clusters,be);
    B.cellnumber = X.cells(2);
    GetFigure('CompareOneShape','parent',X.parent);
    CompareSpikeShape(A,B);
    
function CompareCallback(src, b, cells)
    DATA = GetDataFromFig(src);
    PC.CompareCells(DATA,[],DATA.CellList,cells(1),cells(2));
    
function CombineCells(src, b, cella, cellb)

    DATA = GetDataFromFig(src);
    DATA.CellList = celllist.CombineCells(DATA.CellList, cella, [cella cellb]);
    PC.PlotCellList(DATA);
    PC.CompareCells(DATA,[],DATA.CellList,cella,cella);
    SetData(DATA);

function FindPairs(src,b)
    
    DATA = GetDataFromFig(src);
    C = DATA.CellList;
    mbar = get(get(src,'parent'),'parent');
    cells = unique(C(C>0));
    pairs = [];
    for j = 1:length(cells)
        id = find(C == cells(j));
        [celle{j},cellp{j},cellcl{j}] = ind2sub(size(C),id);
        for k = 1:j-1;
            e = bsxfun(@minus,celle{j},celle{k}');
            p = bsxfun(@minus,cellp{j},cellp{k}');
            id = find(abs(e) < 3 & p ==0);
            if ~isempty(id)
                [a,b] = ind2sub(size(p),id);
                pairs(end+1,:) = [cells(j) cells(k)];
                probes(size(pairs,1)) = cellp{j}(a(1));
            end
        end
    end
    it = findobj(allchild(mbar),'flat','tag','CheckPairs');
    if isempty(it)
        it = uimenu(mbar,'label','Compare','tag','CheckPairs');
    else
        delete(allchild(it));
    end
        
    for j = 1:size(pairs,1)
        uimenu(it,'label',sprintf('P%d %d & %d',probes(j),pairs(j,:)),'callback',{@CompareCallback, pairs(j,:)});
    end
    