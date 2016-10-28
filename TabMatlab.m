function txt = TabMatlab(name, outname)
%TabMatlab(name, @class) adds spaces to align a .m file
%meant for class def files (like @PC/PC.m) to align function names
%
if ~exist(name,'file')
    name = which(name);
end
txt = scanlines(name);

methid = find(CellToMat(strfind(txt,'method')));
endid = find(CellToMat(regexp(txt(methid(1):end),'^\s*end')));
lid  = (1+methid(1)):(methid(1)+endid(1)-2);
eqid = regexp(txt,'=');
%need a factor of 1.1 to make spaces match letters
pos = max(cat(1,eqid{lid}));
spadj = 0.5;
spn = floor(pos * (1+spadj));
good = ones(1,length(txt));
for j = lid(:)'
    l = 1+j-lid(1);
    nsp = 0;
    if isempty(txt{j})
        good(j) = 0;
        sortstr{j} = '';
    elseif isempty(eqid{j})
        id = regexp(txt{j},'\S');
        if ~isempty(id)
            nsp = 2+spn-id(1);
        end
        sortstr{j} = lower(txt{j}(id(1):end));
    else
        id = regexp(txt{j}(1:eqid{j}),'\S');
        sortstr{j} = lower(strtrim(txt{j}(1+eqid{j}:end)));
        nc = length(id);
        nx =  sum(ismember(txt{j}(1:eqid{j}),',[]'));
        id = regexp(txt{j}(1:eqid{j}),'[A-Z]');
        nU = length(id);
        a = eqid{j}(1);
        nsp = spn -a - round(nc * spadj + nU * 0.3 - nx * 0.5);
    end
    if nsp > 0
        addstr = repmat(' ',1,nsp);
        txt{j} = [addstr txt{j}];
    end
end
sortstr{length(txt)} = '';
sortstr = sortstr(find(good));
txt = txt(find(good));
endid = find(CellToMat(regexp(txt(methid(1):end),'^\s*end')));
lid  = (1+methid(1)):(methid(1)+endid(1)-2);
[~,b] = sort(sortstr(lid));
txt(lid) = txt(lid(b));

outname = strrep(name,'.m','.new');
if ~strcmp(outname,name)
    writelines(outname,txt);
end