function DATA = Check(DATA, type, varargin)
%DATA = PCCheck(DATA, type, varargin) makes sure
%status of DATA is useable. type determines what is checked
%'selectprobe'   DATA.selectprobe /selectautoprobe

listtype = PC.GetValue(DATA,'listtype');
CellList = PC.GetValue(DATA,'CellList');
if strcmp(type,'selectprobe')
    if strcmp(listtype, 'autolist')
        if ~isfield(DATA,'selectautoprobe') || size(DATA.selectautoprobe,1) < size(CellList,1) ...
                || size(DATA.selectautoprobe,2) < size(CellList,2) 
            DATA.selectautoprobe(size(CellList,1),size(CellList,2)) = 0;
        end
    else
        if ~isfield(DATA,'selectprobe') || size(DATA.selectprobe,1) < size(CellList,1) ...
                || size(DATA.selectprobe,2) < size(CellList,2) 
            DATA.selectprobe(size(CellList,1),size(CellList,2)) = 0;
        end        
    end
end
        