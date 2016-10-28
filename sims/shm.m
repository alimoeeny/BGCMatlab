function shm()

[t,x] = ode45(@yourfcn,[0 5], [0 1]);
y = x(:,1);
plot(t,y);

function dx = yourfcn(t,x)
w=6.21;
dx(1)=x(2);
dx(2)=-w*sin(x(1));
dx =dx';