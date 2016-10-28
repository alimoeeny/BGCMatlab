function [S, K] = RemoveRedundantRows(S, cmaxcrit,idstr)
%[S, K] = RemoveRedundantRows(S, cmaxcrit, idstr) 
%calculates correlation coefficients. Removes rows from S
%that are correlated with other rows > cmaxcrit

if nargin < 3
    idstr = '';
end
flipped = 0;
plotalong = 0;
if diff(size(S)) > 0
    C = cov(S');
    flipped = 1;
else
    C = cov(S);
end
for j = 1:size(C,1)
    for k = 1:j-1
        ccoeff(j,k) = C(j,k)./(sqrt(C(j,j) * C(k,k)));
        ccoeff(k,j) = ccoeff(j,k);
    end
end

[cmax,b] = max(abs(ccoeff(:)));
K.corrmax(1) = cmax; %largeest at start
rmid = [];
corrs = cmax;
while cmax > cmaxcrit
    [a,b] = ind2sub(size(ccoeff),b);
    if plotalong
        plot(S(:,a),S(:,b),'.');
        title(sprintf('%.2f',ccoeff(a,b)));
    end
    ccoeff(a,b) = 0;
    ccoeff(b,a) = 0;
    rowscores = sum(abs(ccoeff(:,[a b])));
    maxscores = max(abs(ccoeff(:,[a b])));
    if diff(rowscores) > 0 && diff(maxscores) > 0
        rmid(end+1) = b;
    elseif diff(rowscores) < 0 && diff(maxscores) < 0
        rmid(end+1) = a;
    elseif rowscores(1) .* maxscores(1)./(rowscores(2).*maxscores(2)) > 1
%element 1 of one pair is biggest
            rmid(end+1) = a;
    else
            rmid(end+1) = b;
    end
    if rmid(end) == a
        K.keptid(length(rmid)) = b;
    else
        K.keptid(length(rmid)) = a;
    end
    corrs(length(rmid)) = cmax;
    ccoeff(rmid(end),:) = 0;
    ccoeff(:,rmid(end)) = 0;
    [cmax,b] = max(ccoeff(:));
end
K.corrmax = corrs; %largest at end
K.corrmax(end+1) = cmax;
K.rmid = rmid;
K.good = setdiff(1:size(C,1),rmid);
if ~isempty(rmid)
    fprintf('%sEliminating %d redundant Rows of %d:',idstr,length(rmid),length(C));
    for j = 1:length(rmid)
        fprintf(' %d',rmid(j));
    end
    fprintf(' next %.2f\n',cmax);
    S(:,rmid) = [];
end

if flipped
    S = S';
end