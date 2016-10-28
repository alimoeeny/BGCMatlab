function S = Merge(A,B)

if isempty(A);
    S = B;
elseif isempty(B)
    S = A;
else
    S = A;
    if isstruct(S)
    try 
        S(end+1:end+length(B)) = B;
    catch
        nt = length(S);
        f = fields(B);
        for j = 1:length(B)
            for k = 1:length(f)
                S(nt+j).(f{k}) = B(j).(f{k});
            end
        end
    end
    elseif iscell(S)
    end
end
