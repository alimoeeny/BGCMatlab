function result = regressbias(varargin)

shape = 'linear';
nrep = 1000;
noises = 1:5;
amps = 1;
usetan = 0;
std = 3;
nvals = 10;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'amp',3)
        j = j+1;
        amps = varargin{j};
    elseif strncmpi(varargin{j},'gauss',3)
        shape = 'gauss';
    elseif strncmpi(varargin{j},'noise',3)
        j = j+1;
        noises = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        std = varargin{j};
    elseif strncmpi(varargin{j},'nrep',3)
        j = j+1;
        nrep = varargin{j};
    end
    j = j+1;
end
hold off;

if strcmp(shape,'gauss')
    x = -3:3;
    vals = Gauss(1,x);
    vals = vals .* 6./max(vals);
else
    vals = [1:1:nvals]-1;
end

for ai = 1:length(amps)
    vals = amps(ai) .* vals./max(vals);
    rnd = randn(nrep,length(vals)).*std;
    xrnd = randn(nrep,length(vals)).*std;
    for k = 1:length(noises)
        for j = 1:size(rnd,1);
            y = vals + rnd(j,:);
            x = vals + xrnd(j,:)./noises(k);
            [a,b] = fit_bothsubj2error(x,y,noises(k).^2);
            xc = corrcoef(x,y);
            corrs(j) = xc(1,2);
            slopes(j) = b;
            allx(j,:) = x;
            ally(j,:) = y;
        end
        result.corrmean(k,ai) = mean(corrs);
        result.slopem(k,ai) = median(slopes);
        [a,b,c] = Rayleigh(atan((slopes)),'axial');
        result.slopemean(k,ai) = tan(c);
        result.tanmean(k,ai) = tan(mean(atan(slopes)));
        [a,b] = hist(slopes,[0:0.1:4]);
        plot(b,a);
        hold on;
    end
end

result.slopes = slopes;
result.vals = vals;

