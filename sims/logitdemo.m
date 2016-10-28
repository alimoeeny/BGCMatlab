function logitdemo(theta)

mx = 5;
b = [1 1];
[x,y] = meshgrid([-mx:0.01:mx],[-mx:0.01:mx]);
r = x.* b(1).* cos(theta) + y .* b(2).* sin(theta);
hold off;
imagesc(minmax(x(:)),minmax(y(:)),erf(r));
hold on;
asamples = randn(10,2);
asamples(:,2) = asamples(:,2)+0.5; 
asamples(:,1) = asamples(:,1)+0.5; 
plot(asamples(:,1),asamples(:,2),'yo','markerfacecolor','r');
bsamples = randn(10,2);
bsamples(:,2) = bsamples(:,2)-0.5; 
bsamples(:,1) = bsamples(:,1)-0.5; 
plot(bsamples(:,1),bsamples(:,2),'wo','markerfacecolor','b');