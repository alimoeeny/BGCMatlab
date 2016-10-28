function  [score, details] = CalcFitScore(fit)
%[score, details] = CalcFitScore(fit)
%estimate how many isolated cells this found
if ~isfield(fit,'isolation') || isempty(fit.isolation)
    score = NaN;
    details.scores = [NaN NaN NaN];
    return;
end
a = sum(fit.isolation(:,1) > 3);
if size(fit.isolation,2) > 2
    b = sum(fit.isolation(:,3) > 2);
    c = sum(fit.isolation(:,3) > 3 & fit.isolation(:,3) > 2);
    score = a + b + c;
    id = find(fit.isolation(:,1) > 2.5);
    is = (fit.isolation(id,1) -2.5).^0.3;
    score = score + sum(is)./10;
else
    score = a;
    c = 0;
    is = 0;
end

details.scores(1) = score;
details.scores(2) = a+c;
details.scores(3) = sum(is);