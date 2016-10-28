function cell = GetCellNumber(name)
%Get a cell number from a string,Expt or Cluster
%not all implemented yet!
cell = NaN;

if iscell(name)
    for j = 1:length(name)
            cell(j) = GetCellNumber(name{j});
    end
elseif expt.isall(name)
    cell = unique([name.Header.cellnumber]);
    cell = cell(cell > 0);
elseif isfield(name, 'Header') && isfield(name.Header,'cellnumber')
    cell = name.Header.cellnumber;
elseif isfield(name, 'Header') && isfield(name.Header,'loadname')
    cell = NumberFromString(name.Header.loadname, 'cell');
elseif isfield(name, 'Data') && isfield(name.Data,'Header')
    H = name.Data.Header;
    if isfield(H,'cellnumber')
        cell = H.cellnumber;
    elseif isfield(H,'loadname')
        cell = NumberFromString(H.loadname, 'cell');
    else
        cell = NumberFromString(name.name, 'cell');
    end
elseif isfield(name, 'name') 
    cell = NumberFromString(name.name, 'cell');
elseif ischar(name)
    cell = NumberFromString(name, 'cell');
end
if isempty(cell) || sum(~isnan(cell) ==0)
    cell = NaN;
end

