function pc = Afraz(varargin)
%look at effects of slope/bias on % correct at 84/46

signal = -2:0.01:2;
y = erf(signal);
z = erf(signal./1.05);

tsi = 104;
si = [201 - tsi  201 + tsi]
pc(1) = (abs(y(si(1))) + y(si(2)))/2;
pc(2) = (abs(z(si(1))) + z(si(2)))/2;
bias = 0.1;
b = erf(signal+bias);
pc(3) = (abs(b(si(1))) + b(si(2)))/2;
while pc(3) > 0.84
    bias = bias + 0.01;
    b = erf(signal+bias);
    c = erf(signal-bias);
    pc(3) = (abs(b(si(1))) + b(si(2)))/2;
end

GetFigure('SimPsych');
hold off; plot(signal,y); hold on; plot(signal,z);
plot(signal,b);
plot(signal,c);