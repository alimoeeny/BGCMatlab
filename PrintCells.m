function PrintCells(C, varargin)
%PrintCells(C) %prints a cell string array
%PrintCells(C,'tofile',filename)
%(more compact than format compact
%
matchpat = [];
j = 1;
fid = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'match',5)
        j = j+1;
        matchpat = varargin{j};
    elseif strncmpi(varargin{j},'tofile',5)
        j = j+1;
        filename = varargin{j};
        fid = fopen(filename,'w');
        if fid < 1
            fid = 1;
        end
        
    end
    j = j+1;
end

if iscell(C)
    for  j = 1:length(C)
        if isempty(matchpat) | regexp(C{j},matchpat)
            fprintf(fid,'%s\n',C{j});
        end
    end    
elseif ischar(C)
    fprintf(fid,'%s\n',C);
end

if fid > 1
    fclose(fid);
end