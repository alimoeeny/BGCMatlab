function p = anova1u(disptuning)
% ANOVA test whether it is tuned to disparity
% See help for anovan. Although I am only doing one-way anova, I have to use anovan because we may have unequal reps
X = []; group{1} = [];
nd = length([disptuning.x]);
if isfield(disptuning,'counts')
    counts = disptuning.counts;
elseif isfield(disptuning,'rates')
    counts = disptuning.rates;
end
for jj=1:nd
    if isfield(disptuning,'n')
        n = disptuning.n(jj);
    else
        n = length(counts{jj});
    end
    X = [ X sqrt(counts{jj}) ];
    group{1} = [ group{1} [disptuning.x(jj)].*ones(1,n)];
end

[p,table] = anovan(X,group,[],[],[],'off');

