function Both = CatStruct(To, From, varargin)
% Both = CatStruct(To, From) Add struct elements even if fields don't match
%see also CopySFields
% Both = CatStruct(To, From, ids) copies From(j) to To(ids(j));

toid = [];
j = 1;
while j <= length(varargin)
    if isnumeric(varargin{j})
        if length(varargin{j}) == length(From)
            toid = varargin{j};
        end
    end
    j = j+1;
end

if isempty(From)
    Both = To;
    return;
end
ns = length(To);
f = fields(From);
if isempty(toid)
    toid = ns+1:ns+length(From);
end
for j = 1:length(From)
    for k = 1:length(f)
        To(toid(j)).(f{k}) = From(j).(f{k});
    end
end
Both = To;

