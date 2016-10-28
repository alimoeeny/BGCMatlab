function ratio = ACsimple(varargin)
%simple test of how AC ratio with OPNL depends on input strength.

nt = 100000;
sigma = 0.001;
p = 2;
trunc = 0;

x =  randn(nt,1) .* sigma;
nblank = round(nt.*0.9);
%nblank = 0;
if nblank > 0 %make a bunch of frames zero at random
    x(1:nblank) = 0;
end
y = x(randperm(length(x)));
ur = (x + y).^2;
br = (x+x).^2;

%binocular response 
umean(1) = mean(ur);
bmean(1) = mean(br);
umean(2) = mean(ur.^p);
bmean(2) = mean(br.^p);
ratio = umean./(bmean-umean);
xratio = bmean./umean
