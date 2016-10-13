function [disps, dresp, allresps, antiresps, details] = rls(nreps, varargin)
%[disps, dresp, allresps, antiresps, details] = rls(nreps)
%simple simulation of 1D BEM with 1D noise stimulus
%rls(n, 'norm', 'gauss|TI|quad') sets method for monocular normalization
%  monoc resp = X./(1+kN^cpower) where X is RF x Stim, and N is
%  sqrt(E^2+O^2).  Default k = 0.1, cpower = 1.5
% to set these rls(n,'norm','quad',[k cpower]);
showplot = 0;

pix2deg = 0.01;
pix2deg = 0.1;
contrast = [1 1];
normalize = [0 0];
ngain = [0 0 0];
ralfplot = 0;
showhists = 1;
dx = 0;
showrf = 0;
noisetype = 'binary';
binary = 1;
seed = 0;
cpower = 1.5; 
showerr = 0;
getrf = 0;
%ncoeff = 1, nfixed = 0 = perfect contrast gain control, so that a^2+b^1 =
%1 in each eye. 
ncoeff = 0.1;
nfixed = 1;
normalization = 'gauss';
normalization = 'TI';
%Need to make defualt normalization 'none' for allresp to be sensible
normalization = 'none';
npix = 256;
stim = 'rls';
dotw = 0;


sd = 0.4;
rnd = [];
disps = [];
sf = 1;
j = 1;
while j < nargin
    if strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'disps',5)
        j = j+1;
        disps = varargin{j};
    elseif strncmpi(varargin{j},'dotw',4)
        j = j+1;
        dotw = varargin{j};
    elseif strncmpi(varargin{j},'dtplot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'dx',2)
        j = j+1;
        dx = varargin{j}/2;
    elseif strncmpi(varargin{j},'gauss',5)
        noisetype = 'Gauss';
        binary = 0;
    elseif strncmpi(varargin{j},'lc',2)
        j = j+1;
        contrast(2) = varargin{j};
    elseif strncmpi(varargin{j},'getrf',4)
        getrf = 1;
    elseif strncmpi(varargin{j},'norm',4)
        j = j+1;
        normalization = varargin{j};
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            x = varargin{j};
            ncoeff = x(1);
            if length(x) > 1
                cpower = x(2);
            end
        end
    elseif strncmpi(varargin{j},'npix',4)
        j = j+1;
        npix = varargin{j};
    elseif strncmpi(varargin{j},'pix',3)
        j = j+1;
        pix2deg = varargin{j};
    elseif strncmpi(varargin{j},'seed',4)
        j = j+1;
        seed = varargin{j};
    elseif strncmpi(varargin{j},'sqwave',4)
        j = j+1;
        stim = 'sqwave';
    elseif strncmpi(varargin{j},'rc',2)
        j = j+1;
        contrast(1) = varargin{j};
    elseif strncmpi(varargin{j},'rfplot',4)
        showrf = 1;
        showplot = 1;
        if strncmpi(varargin{j},'rfplotonly',8)
            showrf = 2;
        end
    elseif strncmpi(varargin{j},'ngain',4)
        j = j+1;
        ngain = varargin{j};
    elseif strncmpi(varargin{j},'rnd',2)
        j = j+1;
        rnd = varargin{j};
    elseif strncmpi(varargin{j},'sf',2)
        j = j+1;
        sf = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sd = varargin{j};
    end
    j = j+1;
end

if strcmp(normalization,'none')
    ncoeff = 0;
end
Ra = Gabor([sf sd 0 1 -dx],'pix2deg',pix2deg,'npts',npix);
La = Gabor([sf sd 0 1 dx],'pix2deg',pix2deg,'npts',npix);
Rb = Gabor([sf sd pi/2 1 -dx],'pix2deg',pix2deg,'npts',npix);
Lb = Gabor([sf sd pi/2 1 dx],'pix2deg',pix2deg,'npts',npix);
Rc = Gabor([0 sd 0 1 -dx],'pix2deg',pix2deg,'npts',npix);
Lc = Gabor([0 sd 0 1 dx],'pix2deg',pix2deg,'npts',npix);

rfscale = sqrt(sum([Ra.^2 Rb.^2 La.^2 Lb.^2]));
if getrf %show rf
    disps = Ra;
    dresp = rfscale;
    allresps = Rb;
    return;
end
Ra = Ra./rfscale;
Rb = Rb./rfscale;
La = La./rfscale;
Lb = Lb./rfscale;


if showrf
    x = ([1:length(Ra)]-length(Ra)/2).*pix2deg;
    GetFigure('BEM RFs');
    hold off;
    plot(x,Ra,'c');
    hold on;
    plot(x,Rb,'g');
    plot(x,La,'b');
    plot(x,Lb,'r');
    plot(x,Rc,'g:');
    plot(x,Lc,'r:');
    if showrf == 2
        return;
    end
end


normscale(1) = sqrt(1+ ngain(1) * contrast(1).^2);
normscale(2) = sqrt(1+ ngain(2) * contrast(2).^2);
normscale(3) = sqrt(1+ngain(3) * mean(contrast.^2));


%make all disparities multiples of two pixels
nd = 1;
if ~isempty(disps)
    disps = 2.*round(disps./(pix2deg *2));
    maxdisp = max(disps);
else
    maxdisp = 20;
    disps = -(maxdisp):2:(maxdisp);
end
maxdisp = max([1 maxdisp]);

npix = npix + maxdisp*2;
if seed > 0
    rng(seed);
end


if strcmp(stim,'sqwave')
    rnd = BuildSqwave(npix+maxdisp.*2,pi/2);
end
if isempty(rnd)
    if strcmp(noisetype,'Gauss')
        rnd = randn(npix,nreps);
    else
        rnd = rand(npix,nreps) - 0.5;
    end
end

if dotw > 0
    rnd(find(rnd < 0)) = -1;
    rnd(find(rnd > 0)) = 1;
    dotrnd = ceil(rand(1,size(rnd,2)).*dotw);
    for j = 1:dotw:npix
        id = j:j+dotw-1;
%don't make all teh bar transitions at the same locations in 
%every stim.
        for k = 1:size(rnd,2)
            pid = id+dotrnd(k);
            rnd(pid,k) = rnd(pid(1),k);
        end
    end
    rnd = rnd(1:npix,:);
elseif binary
    rnd(find(rnd < 0)) = -1;
    rnd(find(rnd > 0)) = 1;
end
if ralfplot
    GetFigure('PosSums');
    hold off;
    GetFigure('PhaseSums');
    hold off;
    colors = mycolors;
end

if strcmp(normalization,'gauss')
    lnorm = zeros(length(disps),nreps);
    rnorm = zeros(length(disps),nreps);
end
for disp = disps(:)';
    for j = 1:nreps
        iml = contrast(2) .* rnd(maxdisp-disp/2:maxdisp-disp/2+length(Ra)-1,j)'./normscale(2);
        imr = contrast(1) .* rnd(maxdisp+disp/2:maxdisp+disp/2+length(Ra)-1,j)'./normscale(1);
        
        rresp(1,j) = sum(Ra .* imr);
        lresp(1,j) = sum(La .* iml);
        rresp(2,j) = sum(Rb .* imr);
        lresp(2,j) = sum(Lb .* iml);
        if ncoeff > 0
            if strcmp(normalization,'gauss')
                rnorm(nd,j) = sum((Rc(1:end-1).*diff(imr)).^2);
                lnorm(nd,j) = sum((Lc(1:end-1).*diff(iml)).^2);
                ln = lnorm(nd,j).^cpower;
                rn = rnorm(nd,j).^cpower;
            else
                ln = sqrt(sum(lresp(:,j).^2))^cpower;
                rn = sqrt(sum(rresp(:,j).^2))^cpower;
            end
            lnresp(:,j) = lresp(:,j)./(nfixed + ncoeff.* ln);
            rnresp(:,j) = rresp(:,j)./(nfixed + ncoeff.* rn);
        else
        end
    end
    if ncoeff > 0
        allimresps(nd,:) = cat(1,lnresp(:,end), rnresp(:,end));
    end
    
    aresp = (lresp - rresp)./normscale(3);
    resp = (lresp + rresp)./normscale(3);
    qresp = (lresp + flipud(rresp));
    tiresp = (lresp-rresp).^2;
    sqresp = resp .^2;
    dresp(nd) = mean(sum(sqresp,1));
    if strcmp(normalization,'TI')
        nresp = sqresp./(1+tiresp);
    elseif strcmp(normalization,'none')
        nresp = resp.^2;
    else
        nresp = (lnresp + rnresp)./normscale(3);
    end

    nsqresp = nresp .^2;
%sum subunits. Take mean over stimuli    
    ndresp(nd) = mean(sum(nsqresp,1));
    disps(nd) = disp;
    
    
    allresps(nd,:) = sum(sqresp);
    if nargout > 4
        if strcmp(normalization,'none')
            evenresps(nd,:) = resp(1,:);
            oddresps(nd,:) = resp(2,:);
        else
            evenresps(nd,:) = nresp(1,:);
            oddresps(nd,:) = nresp(2,:);
        end
    end
    dsd(nd) = std(allresps(nd,:));
    vmr(nd) = var(allresps(nd,:))./mean(allresps(nd,:));
    antiresps(nd,:) = sum(aresp.^2);
    if ralfplot
        GetFigure('PosSums')
        c = 1+mod(nd-1,size(colors,2))
        plot(rresp(1,:)+lresp(1,:),rresp(2,:)+lresp(2,:),'.','color',colors{c});
        hold on;
        GetFigure('PhaseSums')
        c = 1+mod(nd-1,size(colors,2))
        plot(rresp(1,:)+lresp(2,:),-rresp(2,:)+lresp(1,:),'.','color',colors{c});
        hold on;
    end
    allresp(nd,:) = sum(sqresp);
    allnresp(nd,:) = sum(nsqresp);
    nd = nd+1;
end

disps = disps .* pix2deg;
if nargout > 4
    details.evenresps = evenresps;
    details.oddresps = oddresps;
    details.rf = [La; Ra; Lb; Rb;];
    if ncoeff > 0
        details.lnresp = lnresp;
        details.rnresp = rnresp;
    end
end
if showplot
    GetFigure('BEM Disparity Tuning');
    hold off;
    dscale = mean(dresp([1 end]));
    if showerr
        errorbar(disps,dresp./dscale,dsd./dscale,'color','k');
        hold on;
        plot(disps,dresp./dscale,'ko');
    else
        plot(disps,dresp./dscale,'ko-');
        hold on;
    end
    if ncoeff > 0
        dscale = mean(ndresp([1 end]));
        plot(disps,ndresp./dscale,'ro');
    end
end

details.resp = resp;
details.rnd = rnd;
zid = find(disps == dx); %pref disp
if nreps == 1
    GetFigure('BEM Disparity Tuning');
    hold off;
    plot(disps,ndresp./ndresp(zid));
    hold on;
    plot(disps,dresp./dresp(zid));
    [a,b] = max(ndresp);
elseif showhists
    GetFigure('BEM Histograms');
    hold off;
    did(1) = 1;
    [~,did(2)] = max(dresp);
    [~,did(3)] = min(dresp);
    
    for j = did(:)'
        [y,x] = smhist(allresp(j,:),'noclip','sd',3); plot(x,y);
        hold on;
    end
    GetFigure('Normalised BEM Histograms');
    sd = max(allnresp(:) ./ sqrt(nreps));
    hold off;
    for j = did(:)'
        [y,x] = smhist(allnresp(j,:),'sd',sd,'noclip'); plot(x,y);
        hold on;
    end
end
if ~strcmp(normalization,'none');
     allresps = allnresp;
end
details.allresps = allresps;
details.allnresp = allnresp;

function x = BuildSqwave(npix,dp)

f = 1:2:40;
f = f.*3;
for j = 1:length(f)
    x(j,:) = sin(([1:npix] * 2 .*pi .*f(j)./npix)+dp)./f(j);
end
x = sum(x)';

