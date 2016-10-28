function comp = CompareCellsShape(DATA, CellList, Clusters, cid, varargin) 
%comp = PC.CompareCells(DATA, CellList, Clusters, cid, varargin)     
%Builds a matrix of shape cross correlations
    plottype = '';
    checkmerge = 0;
    if isfield(DATA,'ArrayConfig')
        Array = DATA.ArrayConfig;
    else 
        Array = [];
    end
    j = 1;
    while j <= length(varargin)
        if strncmp(varargin{j},'merged',4);
           checkmerge = 1;
        elseif strcmpi(varargin{j},'noplot')
            plottype = 'none';
        elseif strncmp(varargin{j},'plot',4);
           plottype = 'full';
        end
        j = j+1;
    end
    [ea, pa, ca] = celllist.find(CellList,cid(1));
    [eb, pb, cb] = celllist.find(CellList,cid(2));
    xc = [];
    efficacy = [];
    proberange = PC.GetValue(DATA,'shapecorr');
    
    for j = 1:length(ea)
        for k = 1:length(eb)
            C{1} = PC.GetClusterInfo(Clusters,[ea(j) pa(j) ca(j)],'allexpt');
            C{2} = PC.GetClusterInfo(Clusters,[eb(k) pb(k) cb(k)],'allexpt');
            [xc(j,k), details] = PC.ShapeCorr(C{1},C{2},'delays',5,'drifts',2,Array,'proberange',proberange);
            if xc(j,k) < 0.8
            end
            if ea(j) == eb(k)  && checkmerge == 0%same expt, check efficacy
                e = CalcEfficacy(C{1},C{2});
                efficacy(ea(j)) = max(e(1:2));
            end
        end
    end
    comp.p = prctile(cat(1,pa, pb),50);
    comp.cells = cid;
    comp.ea = ea;
    comp.eb = eb;
    comp.pa = pa;
    comp.pb =pb;
    comp.ca = ca;
    comp.cb = cb;
    comp.xc = xc;
    comp.nonmatch = [];
    comp.goodmatch = [];
    comp.merge = -1;
    comp.efficacy = efficacy;
    eboth = intersect(ea,eb);
    gid = find(xc > 0.95);
    bid = find(xc < 0.8);
    xcm = mean(xc(:));
    if checkmerge
        comp.merge = 1;
        xca = mean(xc,2);
        id = find(xca < 0.9);
        comp.nonmatch = ea(id);
        if ~isempty(id)
            comp.merge = 2;
        end
        if size(xc,1) > 1
            transitions = diag(xc,1);
            id = find(transitions < 0.9);
            comp.badstep = ea(id);
        else
            comp.badstep = [];
        end
%if a whole expt it bad, guarantees that some transitions will be        
        if ~isempty(id) && comp.merge == 0
            comp.merge = 3;
        end
        return;
    end
    if ~isempty(gid)
        did = find(efficacy(eboth) < 0.1);
        sid = find(efficacy(eboth) > 0.2);
        for j = 1:length(eboth)
            
        end
        if ~isempty(did)
            if length(sid) > length(did)
                fprintf('Cells %d and %d xc%.3f could merge %.3f %d/%d sync, %d not\n',cid,xcm,length(sid),length(eboth),length(did));
                comp.merge = 3
            else
                fprintf('Cells %d and %d: xc%.3f. %d/%d efficacy < 0.1, %d > 0.2\n',cid,xcm,length(did),length(eboth),length(sid));
                comp.merge = -1;
                if ~isempty(sid) %need to check wehre it appeas to be the same; 
                    comp.merge = -2;
                end
            end
        elseif ~isempty(bid)
            fprintf('Cells %d and %d could merge %.3f %d good, but %d bad',cid,xcm,length(gid),length(bid));
            [a,b] = ind2sub(size(xc),bid);
            nx = [length(unique(a)) length(unique(b))];
            if nx(1) == 1 && size(xc,1) > 2
                fprintf('E%dP%d cell %d is problem',ea(unique(a)),pa(unique(a)),cid(1));
            end
            if nx(2) == 1 && size(xc,2) > 2
                fprintf('E%dP%d cell %d is problem',eb(unique(b)),pb(unique(b)),cid(2));
            end
            fprintf('\n');
            xca = mean(xc,2);
            id = find(xca < 0.9);
            if ~isempty(id)
                fprintf('Cell%d bad at %s',cid(1),sprintf('%d ',ea(id)));
                comp.nonmatch = [comp.nonmatch ea(id)'];
            end
            id = find(xca > 0.95);
            if ~isempty(id)
                comp.goodmatch = [comp.goodmatch ea(id)'];
            end
            xcb = mean(xc,1);
            id = find(xcb < 0.9);
            if ~isempty(id)
                fprintf('Cell%d bad at %s',cid(2),sprintf('%d ',eb(id)));
                comp.nonmatch = [comp.nonmatch eb(id)'];
            end
            id = find(xcb > 0.95);
            if ~isempty(id)
                comp.goodmatch = [comp.goodmatch eb(id)'];
            end
            if xcm > 0.9 || prctile(xc(:),20) > 0.9
                fprintf('Should Merge\n');
                comp.merge = 2;
                if ~isempty(comp.nonmatch)
                    comp.merge = 4;
                end
            else
                comp.merge = 0;
            end
        else
            fprintf('Cells %d and%d should merge mean xc: %.2f\n',cid,xcm);
            comp.merge = 1;
        end
    else
        comp.merge = 0;
    end
    if comp.merge <= 0 && ~isempty(comp.goodmatch)
        comp.merge = -2;
    end
    if min(size(xc)) == 1
        comp.steps = 0;
    else
        comp.steps = diag(comp.xc,1);
        if min(comp.steps) < 0.9
            [a,b] = min(comp.steps);
            pos(1) = comp.ea(b);
            pos(2,1) = comp.eb(b+1);
            pos(1,2) = comp.pa(b);
            pos(2,2) = comp.pb(b+1);
            pos(1,3) = comp.ca(b);
            pos(2,3) = comp.cb(b+1);
            PC.CheckAutoList(DATA, 'square',pos);
        end
    end
    PlotComparison(DATA,comp,'full');
    drawnow; 

    
    function F = PlotComparison(DATA, comp,type)
    if strcmp(type,'none')
        F = 0;
        return;
    end
    stepstr = '';
    if length(comp) > 1
        for j = 1:length(comp)
            if ~strcmp(type,'none')
            F = PlotComparison(DATA,comp(j),type);
            end
            if strcmp(type,'fullwait')
                if j == 1
                    setappdata(F,'waitforbutton',1);
                end
                AddNextButton(F);
                a = getappdata(F,'waitforbutton');
                if a
                    uiwait(gcf);
                else
                    fprintf('Got Jump\n');
                    type = 'none';
                    break;
                end
            end
        end
        return;
    end
    F = GetFigure('CompareShapeMatrix','parent',DATA.toplevel);
    hold off; 
    h = imagesc(comp.xc);
    minstep = min(diag(comp.xc,1));
    stepstr = sprintf('Worst Step %.3f',minstep);
    set(gca,'ytick',1:length(comp.ea),'yticklabel',num2str(comp.ea));
    set(gca,'xtick',1:length(comp.eb),'xticklabel',num2str(comp.eb));
    caxis([max([min(comp.xc(:)) 0.5]) 1]);
    colorbar;
    xlabel(sprintf('Cell%dP%d',comp.cells(1),prctile(comp.pa,50)));
    ylabel(sprintf('Cell%dP%d',comp.cells(2),prctile(comp.pb,50)));
    if strncmp(type,'full',4)
        setappdata(gcf, 'CompData',comp);
        set(h,'buttondownfcn',{@HitCompIm});
    end
    if comp.merge > 0
        if isempty(comp.nonmatch)
            title(sprintf('Merge %s',stepstr));
        else
            title(['Merge but check E' sprintf('%d ',comp.nonmatch) stepstr]);
        end
    else
        if isempty(comp.goodmatch)
            title(sprintf('NoMerge %s',stepstr));        
        else
            title(['NoMerge But Check'  sprintf('%d ',comp.goodmatch) stepstr]);        
        end            
    end
        
    
    
   
    function HitCompIm(src,b)

    DATA = GetDataFromFig(src);
    Clusters = PC.GetValue(DATA,'Clusters');
    xy = get(gca,'currentpoint');
    comp = getappdata(gcf,'CompData');
    a = round(xy(1,2));
    b = round(xy(1,1));
    A = PC.GetClusterInfo(Clusters,[comp.ea(a) comp.pa(a) comp.ca(a)]);
    B = PC.GetClusterInfo(Clusters,[comp.eb(b) comp.pb(b) comp.cb(b)]);
    if comp.ea(a) == comp.eb(b) %same expt - look at crosscorr
        e = CalcEfficacy(A.times,B.times);
        xc = xcorrtimes(A.times,B.times);
        GetFigure('xcorr','parent',DATA.toplevel);
        plot(xc);
        title(sprintf('E%d',comp.ea(a)));
    else
        GetFigure('CompareShape','parent',DATA.toplevel);
        PC.CompareSpikeShape({A B});
        [xc, details] = PC.ShapeCorr(A,B, DATA.ArrayConfig,'shifts',5,'drifts',1);
        title(sprintf('Best xc %.3f',xc));
    end
    pos = [comp.ea(a) comp.pa(a) comp.ca(a); comp.eb(b) comp.pb(b) comp.cb(b)];
    PC.CheckAutoList(DATA, 'square',pos);
    PC.CompareQuickSpks(DATA,pos);
    
            