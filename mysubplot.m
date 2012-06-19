function mysubplot(nr,nc,n, varargin)
wscale = 1;
hscale = 1;
vsep = 0.1;
hsep = 0.1;

j = 1;
while j < length(varargin)
    if strncmpi(varargin{j},'width',5)
        j = j+1;
        wscale = varargin{j};
    elseif strncmpi(varargin{j},'tight',5)
        vsep = 0.02;
        hsep = 0.02;
    end
    j = j+1;
end
dw = wscale * hsep./nc;
dh = hscale * vsep./nr;
w = wscale * (1-hsep)./nc;
h = hscale * (1-vsep)./nr;
for j = 1:length(n)
y(j) = (ceil(n(j)./nc)-1) .* (h+dh);
x(j) = (mod(n(j)-1,nc)) .* (w+dw);
end
w = w + max(x)-min(x);
x = min(x)+hsep/(2*nc);
h = h + max(y)-min(y);
y = min(y)+vsep/(2*nr);
subplot('Position',[x 1-y-h w h]);
set(gca,'xticklabel',[],'yticklabel',[]);
