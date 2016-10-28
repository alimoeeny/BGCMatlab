function Fourier(type,varargin)

if nargin == 0
    type = '2dgap';
end
f = 4 .* pi./256;
tf = 16 .* pi./256;

if strcmp(type,'2dgap')
    [xi,yi] = meshgrid(1:256,1:256);
    im = sin(xi .* f + yi.*tf);
    im2 = sin(xi .* f + yi.*tf*4);
    imagesc(fftshift(abs(fft2(im))));
    for j = 2:2:256
        im(j,:) = 0;
    end
    imagesc(fftshift(abs(fft2(im))));
end