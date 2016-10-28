function temporalft(varargin)
%check FT of 150Hz, 75Hz fake, framerates.


t = zeros(1,3000);
xt = t;
rnd = rand(1,300);
t(1:10:3000) = rnd;
xt(1:10:3000) = rnd;
xt(6:10:3000) = rnd;
hold off;

frq =[1:length(t)] * 0.5;
plot(frq,abs(fft(t)));
hold on;
plot(frq,abs(fft(xt)));
set(gca,'xlim',[0 200])
set(gca,'ylim',[0 20])

