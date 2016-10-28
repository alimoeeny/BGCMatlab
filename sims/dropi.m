function dropi(varargin)
%look at dropi estimate from just mean/sd


dropi = [];
x = randn(1000,1);
triggers = 0:0.1:4;

for loop = 1:10;
x = randn(1000,1);
for j = 1:length(triggers(:))
    t = triggers(j);
    y = x + t;
    y = y(y>0);
    dropi(j,loop) = mean(y)./std(y);
end
end

GetFigure('dropi test');
hold off;
plot(triggers,dropi);
truedrop = (dropi.^3 - 2.3).^0.333;
truedrop(dropi.^3 < 2.3) = 0;
fity = (2.3 + triggers.^3).^0.33333;
hold on;
plot(triggers,truedrop);
plot(triggers,fity,'r');
refline(1);