function dup = FixDoubleCells(DATA, varargin)
    CellList = PC.GetValue(DATA,'CellList');
    Clusters = PC.GetValue(DATA,'Clusters');
    
    for e = 1:size(CellList,1)
        [a,b] = Counts(CellList(e,:,:));
        dup = find(a>1 & b > 0)
    end