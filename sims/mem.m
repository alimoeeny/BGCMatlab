function [ER, OR] = mem(varargin)
%simpe motion energy model

drift = 0;
angle = -pi./10;

rf = Gabor([0.1 5 -pi/4 1 0 0 angle 5]);
orf = Gabor([0.1 5 pi./4 1 0 0 angle 5]);

%imagesc(rf);
f = 5./256;
if drift == 0
grating = sin(2.*pi .* f.*[1:256]);
bgrating = sin((2.*pi .* f.*[1:256])+pi);
stim(1:128,:) = repmat(grating,128,1);
stim(129:256,:) = repmat(bgrating,128,1);
end
%imagesc(stim+rf);
for j = 1:256
    ts = circshift(stim,j,1);
    ER(j) = sum(ts(:) .* rf(:)).^2;
    OR(j) = sum(ts(:) .* orf(:)).^2;
end
hold off;
plot(ER);
hold on;
plot(OR);