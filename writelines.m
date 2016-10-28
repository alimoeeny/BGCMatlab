function writelines(name,cellstr, varargin)
%writelines(namn, cellstr) creates an ascii file name
%by printing each element of cellstr onto a different line
%see also scanlines
silent = 0;
mkmatrix = 0;
method = 0;
j = 1;
argon = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'bufsize',4)
        argon = {argon{:} varargin{j} varargin{j+1}};
        j = j+2;
    elseif strncmpi(varargin{j},'matrix',4)
        mkmatrix = 1;
    elseif strncmpi(varargin{j},'fget',4) %use fgets instead of textread
        method = 1;
    elseif strncmpi(varargin{j},'silent',4)
        silent = 1;
    elseif ischar(varargin{j})
    end
    j = j+1;
end

if ~ischar(name)
    if silent == 0
        mycprintf('errors','writelines: Input argument names is not a char\n',name);
    end
    return;
end
fid = fopen(name,'w');
if fid < 0
    mycprintf('errors','writelines: %s exists but fopen fails\n',name);
    return;
end    

for j = 1:length(cellstr)
    fprintf(fid,'%s\n',cellstr{j});
end
fclose(fid);