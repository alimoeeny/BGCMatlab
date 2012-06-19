function handles = DensityPlot(x,y,varargin)

yrange = minmax(y);
xrange = minmax(x);
setydir = 0;
nbins = [50 50];
sx = 0; sy = 0;
handles = [];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'aaa',3)
    elseif strncmpi(varargin{j},'exact',4)
        j = j+1;
        nbins(1) = diff(xrange)+1;
        nbins(2) = diff(yrange)+1;
    elseif strncmpi(varargin{j},'nbins',4)
        j = j+1;
        nbins = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sx = varargin{j}(1);
        sy = varargin{j}(2);
    elseif strncmpi(varargin{j},'ynormal',5)
        setydir = 1;
    elseif strncmpi(varargin{j},'xrange',5)
        j = j+1;
        xrange = varargin{j};
    elseif strncmpi(varargin{j},'yrange',5)
        j = j+1;
        yrange = varargin{j};
    end
    j = j+1;
end
[gx,gy] = meshgrid(-10:10,-10:10);
if sx == 0
    sx= std(x)/10;
end
if sy == 0 
    sy= std(y)/10;
end
    G = exp(-(gx).^2/sx - (gy).^2/sy);
    G = G./sum(G(:));
    z = zeros(nbins(1)+1,nbins(2)+1);
% ignore spikes where both are set to 0
    xtra(1) = diff(xrange)./(nbins(1)-1);
    xtra(2) = diff(yrange)./(nbins(2)-1);
    xi = 1+floor(nbins(1) * (x-xrange(1))/diff(xrange));
    yi = 1+floor(nbins(2) * (y-yrange(1))/diff(yrange));
    for j = 1:length(x)
        z(yi(j),xi(j)) = z(yi(j),xi(j))+1;
    end
    if sx+sy > 0
    z = conv2(z,G,'same');
    end
    handles(1) = imagesc([xrange(1) xrange(2)+xtra(1)],[yrange(1) yrange(2)+xtra(2)],z);
    if setydir == 1
    set(gca,'Ydir','normal');
    end
