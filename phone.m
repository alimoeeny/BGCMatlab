function phone(name, varargin)

phonebook = '/b/bgc/Win/Outlook/contacts.csv';


txt = scanlines(phonebook);
n = 1;
skip = [];
for j = 1:length(txt)
    if isempty(txt{j})
        skip = [skip j];
    else
        while txt{j}(end) ~= ','
            txt{j} = [txt{j} ',' txt{j+n}];
            skip = [skip j+n];
            n = n+1;
        end
    end
    n = 1;
end
txt = txt(setdiff(1:length(txt),skip));
if strcmp(name,'social')
    id = find(CellToMat(regexpi(txt,'Bruce,,Cumming')));
    a = split(txt{id},',');
    b = [a{61}(1:3) a{15}(18:20) a{15}(21:end)];
    fprintf('%s\n',b);
    return;
end
id = find(CellToMat(regexpi(txt,name)));
f = split(txt{1},',');
secondf = [15:18 86:88]; %fields to put at end
for j = 1:length(id)
    s = split(txt{id(j)},',');
    for k = 1:length(s) 
        if ~isempty(s{k}) && ~ismember(k,secondf)
            fprintf('%s ',s{k});
        end
    end
    for k = 1:length(s) 
        if ~isempty(s{k})  && ismember(k,secondf)
            fprintf('%s ',s{k});
        end
    end
    fprintf('\n');
end