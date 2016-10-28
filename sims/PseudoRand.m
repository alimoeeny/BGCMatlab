function PsuedoRand(nt)

nloops = 1000;
nt = round(nt/2);
x = cat(2,ones(1,nt),ones(1,nt)*-1);
for j = 1:nloops
    y = x(randperm(length(x)));
    p(j) = sum(diff(y) ~=0)./length(diff(x));
end
prev = mean(p)