function to = CopyFields(to, from, varargin)

j = 1;
while j <= length(varargin)
    j = j+1;
end

if isempty(from)
    return;
end
f = fields(from);
for j = 1:length(f)
    to.(f{j}) = from.(f{j});
end