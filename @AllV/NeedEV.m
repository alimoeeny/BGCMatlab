function need = NeedEV(C)
%determine if EigenVectors are needed for this classification

need = 0;
if C.space(1) == 1
    need = need+1;
end
for j = 1:length(C.next)
    if isfield(C.next{j},'space') && C.next{j}.space(1) == 1
        need = need+1;
    end
end