function S = CheckStruct(S, varargin)
%S = CheckStruct(S, varargin) find handles in struct
%

name = '';
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'name')
        j = j+1;
        name = varargin{j};
    end
    j = j+1;
end

if iscell(S)
    for j = 1:size(S)
        S{j} = CheckStruct(S{j},varargin{:});
    end
elseif isstruct(S)
    f = fields(S);
    for j = 1:length(f)
        if isstruct(S(1).(f{j})) || iscell(S(1).(f{j}))
            CheckStruct(S(1).(f{j}),'name',[name '.' f{j}]);
        else
            c = class(S(1).(f{j}));
            if strncmp(c,'matlab.graphics',15)
                fprintf('%s.%s %s\n',name,f{j},c);
            end
        end
    end
end