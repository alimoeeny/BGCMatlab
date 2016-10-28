function emCorr(rate, var)
%simulate effecs of EM on corr

for j = 1:1000;
    x = randn(1,1000).* sqrt(rate);
    y = randn(1,1000).* sqrt(rate);
    z = randn(1,1000).* sqrt(var);
    xc = corrcoef(x+z,y+z);
    c(j) = xc(1,2);
end
fprintf('Corr %.4f, SD %.4f',mean(c),std(c));
