function duplist = FixDoubleCells(DATA, varargin)
%PC.FixDoubleCells(DATA, varargin) Find locations with a cell set twice
%for a single expt and resolve
DATA = GetDataFromFig(DATA);
strargs = cell2cellstr(varargin);

    CellList = PC.GetValue(DATA,'CellList');
    Clusters = PC.GetValue(DATA,'Clusters');
    duplist = [];
    for e = 1:size(CellList,1)
        [a,b] = Counts(squeeze(CellList(e,:,:)));
       dup = find(a>1 & b  > 0);
        for j = 1:length(dup)
            duplist(end+1,:) = [e b(dup(j))];
        end
    end
    
    
    
    
    for j = 1:length(duplist)
        e = duplist(j,1);
        celli = duplist(j,2);
        [p,cl] = find(squeeze(CellList(e,:,:)) == celli);
        [alle, allp, allc] = celllist.find(CellList,celli);
        C = PC.GetClusterInfo(DATA,[e p(1) cl(1); e p(2) cl(2)]);
        xc = PC.ShapeCorr(C, DATA.ArrayConfig);
        if xc > 0.98 %Cant use shape
            
        else
            id = find(abs(alle-e) < 5 & abs(alle-e) > 0); %nearby expts
            ds = sort(unique(abs(alle(id)-e)));
            wts = [5 4 3 2 1 0];
            if ds(2) == 2%have 1 and 2
                wts(3:end) = 0;
            end
            nearC = PC.GetClusterInfo(DATA,[alle(id) allp(id) allc(id)]);            
            for k = 1:length(nearC)
                xcs(1,k) = PC.ShapeCorr(C{1},nearC{k},DATA.ArrayConfig);
                xcs(2,k) = PC.ShapeCorr(C{2},nearC{k},DATA.ArrayConfig);
                w(1:2,k) = wts(round(abs(alle(id(k))-e)));
                nspks(k) = nearC{k}.ncut;
                nspkdiff(1,k) = abs(nearC{k}.rate - C{1}.rate);
                nspkdiff(2,k) = abs(nearC{k}.rate - C{2}.rate);
            end
            score = mean(xcs .*w,2)./mean(w,2);
            nsscore = mean(nspkdiff .*w,2)./mean(w,2);
        end
        fprintf('E%dCell%d xc %.2f: %.2f,%.2f nspkdiff %.1f,%.1f',e,celli,xc,score,nsscore);            
        if diff(score) > 0 && diff(nsscore) < 0 %second is best
            duplist(j,3:4) = [p(1) cl(1)]; %clear first
        elseif diff(score) < 0 && diff(nsscore) < 0
            duplist(j,3:4) = [p(2) cl(2)];
        else
            duplist(j,3:4) = [0 0];
            if sum(strcmp('interactive',strargs))
                if diff(p) == 0
                    [a,b] = min(abs(alle(id)-e)); %nearest expt
                    PC.CompareQuickSpks(DATA,[e p(1) cl(1); alle(id(b)) allp(id(b)) allc(id(b))]);
                else
                    PC.CompareQuickSpks(DATA,[e p(1) cl(1); e p(2) cl(2)]);
                end
            end
        end
        fprintf('Replacing %d.%d\n',duplist(j,3:4));
    end
        