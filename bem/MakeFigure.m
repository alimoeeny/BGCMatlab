function varargout = MakeFigure(varargin)
%Make Figures for Disparity Gradient/Slant Simulatinos.
%MakeFigure('DispXSF') Makes Disparity X SF image plot
results = [];
resps = [];
tags{1} = 'DTsep';
fig = 1;
fsz = 12; 
seed = 83;
rebuild = 0;
j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        results = varargin{j};
    elseif strcmp(varargin{j},'rebuild')
        rebuild = 1;
    elseif strcmp(varargin{j},'seed')
        j = j+1;
        seed = varargin{j};
    elseif ischar(varargin{j})
        fig = varargin{j};
    end
    j = j+1;
end

if fig == 2
    load dorisfstim.mat;
    subplot(2,2,1);
    fillpmesh(res.disps,res.dgs,res.resps{end},'plot');
    subplot(2,2,2);
    imagesc(res.sresps{end});

elseif fig == 1
%
% results{5} needs to be dg x dx, for dsf and dor, as in dxdor.mat
% results{4} needs to be dori x ori, for a fixed slant as in oridori.mat
% results{3} needs to be dsf x dori, for a fixed slant, and is oridsf.mat
%
    if isempty(results)
        load dorisims.mat
    end
    
    subplot(3,1,1);
    hold off;
    fillpmesh(results{5}.dgs, results{5}.disps,results{5}.resps{1},'plot');
    axis('square');
    %shading('interp');
    subplot(3,1,2);
    fillpmesh(results{5}.dgs, results{5}.disps,fliplr(results{5}.sresps{1}),'plot');
    %shading('interp');
    axis('square');
    ylabel('disparity (degrees)','fontsiz',fsz);
    xlabel('disparity gradient (degrees/degree)','fontsiz',fsz);
    
    subplot(3,1,3);
    
    odo= PlotSims(results{4});
    oresp = [0 max(odo.sig')];
    sfo = PlotSims(results{3});
    subplot(3,1,3);
    hold off;
    oris = [0 results{4}.oris];
    plot(oris .* 180/pi, oresp,'linewidth',2);
    sresp = [0 max(sfo.sig' .* -1)];
    hold on;
    oris = [0 results{3}.oris];
    plot(oris.*180/pi, sresp,'r','linewidth',2);
    set(gca,'xlim',[0 90],'xtick',[0 45 90],'ytick',[]);
    axis('square');
    xlabel('RF orientation (degrees)','fontsiz',fsz);
    ylabel('Response','fontsiz',fsz);
elseif strcmp(fig,'DispxSFseeds')
    F = GetFigure('SFxDisp');
    pos = get(F,'position');
%Good seeds with 1024 pixels at 0.005 deg/pix are 99,83 67    
    for j = 1:100
        im = MakeFigure('DispxSF','seed',j);
        newf = GetFigure(sprintf('Seed %d',j),'parent',F,'front');
        set(newf,'position',pos);
        imagesc(im');
        title(sprintf('Seed %d',j));
        drawnow;
    end
elseif strcmp(fig,'DispxSF')  
    F = GetFigure('SFxDisp');
    isf = 22;
    if isappdata(F,'PlotData') && rebuild == 0
        xData = getappdata(F,'PlotData');
        im = xData.im;
    else
        sfs = exp(log(0.25):log(1.1):log(16)) .*1.081;
        for j = 1:length(sfs)
            [x, im(:,j)] = rls(1,'sf',sfs(j),'sd', 0.5./sfs(j),'seed',seed,'npix',1024,'pix',0.005,'disps',[-2:0.01:2],'norm','none');
            rfs(:,j) = rls(1,'sf',sfs(j),'sd', 0.5./sfs(j),'npix',1024,'getrf','pix',0.005);
        end
        xData.rfs = rfs;
        xData.im = im;
        xData.sfs = sfs;
        xData.disps = x;
        xData.seed = seed;
        xData.isf = isf;
        setappdata(F,'PlotData',xData);
    end
    subplot(2,1,2);
    hold off;
    imagesc(xData.disps,xData.sfs,im');
    for j = 1:size(im,2)
        im(:,j) = im(:,j)./max(im(:,j));
    end
    imagesc(xData.disps,xData.sfs,im');
    
    nomv = log(min(xData.sfs)) + (log(xData.sfs(isf))-log(min(xData.sfs)))./(log(max(xData.sfs)) - log(min(xData.sfs))); 
    yl = get(gca,'ylim');
    nomv = yl(1) + diff(yl) .* isf./length(xData.sfs) + diff(yl)./(2.* length(xData.sfs));
    Test.SetLogTics(gca,xData.sfs,'Y','tics',[0.25 0.5 1 2 4 8]);
    nomv = Test.SetLogTics(gca,xData.sfs,'Y','value',xData.sfs(isf));

  
    colormap('gray');
    hold on;
    plot(get(gca,'xlim'),[nomv nomv],'r-');
    varargout{1} = im;
    xlabel('RF Position Disparity (degrees)');
    ylabel('Filter SF (cpd)');
    title('Energy model responses across spatial scales');
%    title(sprintf('seed %d\n',seed));
    subplot(2,1,1);
    hold off;
    plot(xData.disps,im(:,22),'r-','linewidth',2); %good for seed 83
    y = sum(im');
    hold on;
    plot(xData.disps,y./max(y),'b-','linewidth',2);
    set(gca,'xtick',[]);
    ylabel('Response');
    title('Energy model responses at one spatial scale');
    GetFigure('RFs used');
    imagesc(xData.rfs);
elseif strcmp(fig,'Sanada')
    F = GetFigure(fig);
    rfa = Gabor([0.3 3 0 1 0 0 0 4]);
    rfb = Gabor([0.3 3 0 1 0 0 0 2 5]);
    rfc = Gabor([0.3 3 0 1 0 0 0 2 -5]);

    subplot(2,3,1);
    imagesc(rfa);
    colormap('gray');
    set(gca,'xtick',[],'ytick',[]);
    axis('image');
    axis('square');

    subplot(2,3,2);
    imagesc(rfa);
    colormap('gray');
    set(gca,'xtick',[],'ytick',[]);
    axis('image');
    axis('square');

    subplot(2,3,3);
    imagesc(rfb+rfc);
    colormap('gray');
    set(gca,'xtick',[],'ytick',[]);
    xv = [0:0.01:1];
    [x,y] = meshgrid(xv,xv);
    R =exp(-(abs((x-y)).^2)./0.3.^2) .* exp(-(abs((x+y)-1).^2)./0.7.^2);
    axis('image');
    axis('square');

    subplot(2,3,6);
    imagesc(R);
    set(gca,'xtick',[],'ytick',[]);
    set(gca,'ydir','normal');
    axis('image');
    axis('square');
    subplot(2,3,4);
    R =exp(-(abs((x-y)).^2)./0.4.^2) .* exp(-(abs((x+y)-1).^2)./0.4.^2);
    imagesc(R);
    set(gca,'xtick',[],'ytick',[]);
    set(gca,'ydir','normal');
    axis('image');
elseif strcmp(fig,'BEMvar')
    F = GetFigure(fig);
    if isappdata(F,'PlotData') && rebuild == 0
        X = getappdata(F,'PlotData');
    else
        dxs = [0 0.125 0.25 0.5 1];
        [dx,y,X.bem,~, details] = rls(10000,'npix',1024,'pix',0.005,'sf',2,'sd',0.5./2,'dotw',50,'disps',dxs); 
        X.dxtuning(1,:) = y;
        X.evenrespA = details.evenresps;
        [dx,y,X.bem,~, details] = rls(10000,'npix',1024,'pix',0.005,'sf',2,'sd',0.5./2,'disps',dxs); 
        X.evenrespB = details.evenresps;
        X.oddrespB = details.oddresps;
        X.dxtuning(2,:) = y;
        
        [a,b] = max(y);
        [c,d] = min(y);
        X.pref = b;
        X.null = d;
        X.dx = dx;
        
        [a,b] = max(y);
%originally was [1 1]        
        [dx,y,X.nbem,~,details] = rls(10000,'norm','quad',[4 1],'npix',1024,'pix',0.005,'sf',2,'sd',0.5./2,'disps',dxs);
        F = GetFigure(fig);
        X.dxtuning(3,:) = y;
        X.evenrespN = details.evenresps;
        X.oddrespN = details.oddresps;
        setappdata(F,'PlotData',X);
    end
    args = {'linewidth' 2};
    roc(1) = CalcCP(X.bem(X.pref,:),X.bem(X.null,:));
    roc(2) = CalcCP(X.nbem(X.pref,:),X.nbem(X.null,:));
    nrow = 3;
    ncol = 2;
    
    subplot(nrow,ncol,1);
    bins = [0.01:0.01:20];
    hold off;
    [x,y] = smhist(X.bem(X.pref,:),'xval',bins,'sd',0.01,'box');
    x = smooth(x,10);
    x = smooth(x,10);
    xx = x;
    yy = y;
    %10x becuase there are 10 samples per boxcar width
    nc = 1;
    counts(nc,1) = sum(x);
    x = 10 .* x./sum(x);
    counts(nc,2) = sum(x);
    nullx = x;
    area(1) = sum(x);
    plot(y,x,args{:});
    hold on;
    [x,y] = smhist(X.bem(X.null,:),'xval',bins,'sd',0.01,'box');
    nc = nc+1;
    counts(nc,1) = sum(x);
    x = 10 .* x./sum(x);
    x = smooth(x,10);
    x = smooth(x,10);
    counts(nc,2) = sum(x);
    plot(y,x,args{:});
    area(2) = sum(x);
    hold on;
    set(gca,'ylim',[0 0.4],'xtick',[],'xlim',[0 4]);

    
    subplot(nrow,ncol,3);
    hold off; 
    plot(y,nullx,args{:});
    hold on;
    plot(y,x,args{:});
    area(2) = sum(x);
    set(gca,'ylim',[0 0.1],'xtick',[],'xlim',[0 4]);

    a = max([prctile(X.nbem(X.pref,:),90) prctile(X.nbem(X.null,:),90)]);
    a = a.*1.1;
%Plot histogram for Normalized responses    
    subplot(nrow,ncol,5);
    bins = a./100:a./200:a;
    sd = a./100;
    hold off;
    [x,y] = smhist(X.nbem(X.pref,:),'xval',bins,'sd',sd,'box');
    x = smooth(x,10);
    nc = nc+1;
    counts(nc,1) = sum(x);
    x = 4 .* x./sum(x);
    counts(nc,2) = sum(x);

    plot(y,x,args{:});
    area(3) = sum(x);
    hold on;
    [x,y] = smhist(X.nbem(X.null,:),'xval',bins,'sd',sd,'box');
    nc = nc+1;
    counts(nc,1) = sum(x);
    x = 4 .* x./sum(x);
    counts(nc,2) = sum(x);
    x = smooth(x,10);

    plot(y,x,args{:});
    hold on;   
    area(4) = sum(x);
    set(gca,'ylim',[0 0.4],'xtick',[],'xlim',[0 max(y)]);

    xlabel('Response');
    ylabel('Relative Frequency');

    
    
%monocular responses    
    subplot(nrow,ncol,2);
    lw = 2;
    bins = [-30:0.01:30];
    aresps = X.evenrespA(X.pref,:) .* std(X.evenrespB(X.pref,:))./std(X.evenrespA(X.pref,:));

    [a,b] = smhist(aresps,'sd',0.1,'box','xval',bins);
    nc = nc+1;
    counts(nc,1) = sum(a);
    a = 10 .* a./sum(a);
    counts(nc,2) = sum(a);
    hold off;
    plot(b,a,'r','linewidth',lw);

    bresps = X.evenrespB(X.pref,:);
    [a,b] = smhist(bresps,'sd',0.1,'box','xval',bins);
    nc = nc+1;
    counts(nc,1) = sum(a);
    a = 10 .* a./sum(a);
    a = smooth(a,10);
    counts(nc,2) = sum(a);
    hold on;
    plot(b,a,'k','linewidth',lw);
    set(gca,'xlim',[-3 3],'xtick',[]);

    nresps = X.evenrespN(X.pref,:) .* std(X.evenrespB(X.pref,:))./std(X.evenrespN(X.pref,:));
    oresps = X.oddrespN(X.pref,:) .* std(X.evenrespB(X.pref,:))./std(X.oddrespN(X.pref,:));
    nresps = nresps .* 1.1;
    [a,b] = smhist(nresps,'sd',0.1,'box','xval',bins);
    nc = nc+1;
    counts(nc,1) = sum(a);
    a = 10 .* a./sum(a);
    a = smooth(a,10);
    a = smooth(a,10);
    counts(nc,2) = sum(a);
    hold on;
    plot(b,a,'b','linewidth',lw);

    
    subplot(nrow,ncol,4);
    hold off;
    bins = [0.01:0.01:20];
    [a,b] = smhist(aresps.^2,'sd',0.01,'box','xval',bins);   
    nc = nc+1;
    counts(nc,1) = sum(a);
    a = 10 .* a./sum(a);
    counts(nc,2) = sum(a);
    a = smooth(a,10);
    a = smooth(a,10);

    plot(b,a,'r','linewidth',lw);
    hold on;
    [c,d] = smhist(bresps.^2,'sd',0.1,'box');
    [a,b] = smhist(bresps.^2,'sd',0.01,'box','xval',bins);
    rnd = ceil(rand(2,length(bresps)).*length(bresps));
    nc = nc+1;
    counts(nc,1) = sum(a);
    a = 10 .* a./sum(a);
    a = smooth(a,10);
    a = smooth(a,10);
    counts(nc,2) = sum(a);

    plot(b,a,'k','linewidth',lw);
    [a,b] = smhist(nresps.^2,'sd',0.1,'box','xval',bins);
    a = 10 .* a./sum(a);
    plot(b,a,'b','linewidth',lw);
    set(gca,'xlim',[0 20],'ylim',[0 0.2],'xtick',[],'xlim',[0 4]);
    
elseif strcmp(fig,'BEMdiagram')
    F = GetFigure(fig);
    clf
    if isappdata(F,'PlotData');
        X = getappdata(F,'PlotData');
        rfs(:,1) = X.rfs(:,X.isf);
    else
        sfs = exp(log(0.25):log(1.1):log(16));
        sf = 2.0;
        pixsz = 0.005;
        [erfs(:,1),b,orfs(:,1)] = rls(1,'sf',sf,'sd', 0.5./sf,'npix',1024,'getrf','pix',pixsz);
        sf=1;
        [erfs(:,2),b,orfs(:,2)] = rls(1,'sf',sf,'sd', 0.4./sf,'npix',1024,'getrf','pix',pixsz);
    end
    x = [1:size(erfs,1)] .* pixsz;
    x = x-mean(x);
    pid = find(abs(x) < 0.8);
    nl = x.^2;
    nl(x<0)=0;
  
    yp = [0.02 0.24 0.56 0.78];
    xp = [0.05 0.3 0.5];
    hs = [0.2 0.1 0.2];
    ws = hs;
    subplot('Position',[0.05 yp(1) 0.2 0.2]);
    h = plot(x(pid),erfs(pid,1),'k');
    set(h,'linewidth',2);
    axis('tight');
    set(gca,'xtick',[],'ytick',[],'color','none');
    
    dy = 0;
    annotation('arrow',[0.25 xp(2)],[0.11 0.23]+dy,'units','normalized');

    subplot('Position',[0.05 yp(2) 0.2 0.2]);
    h=plot(x(pid),erfs(pid,1),'k');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');
    annotation('arrow',[0.25 xp(2)],[0.33 0.23]+dy,'units','normalized');
    
    subplot('Position',[xp(2) dy+0.23-hs(2)/2 ws(2) hs(2)]);
    h = plot(nl,'k');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');
    
    subplot('Position',[0.05 yp(3) 0.2 0.2]);
    h=plot(x(pid),-erfs(pid,1),'k');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');

    subplot('Position',[0.05 yp(4) 0.2 0.2]);
    h=plot(x(pid),-erfs(pid,1),'k');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');

    dy = 0.6;
    annotation('arrow',[0.25 xp(2)],[0.33 0.23]+dy,'units','normalized');
    annotation('arrow',[0.25 xp(2)],[0.11 0.23]+dy,'units','normalized');
    subplot('Position',[xp(2) dy+0.23-hs(2)/2 ws(2) hs(2)]);
    h = plot(nl,'k');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');

    dy=0;
    subplot('Position',[xp(3) dy+(yp(2)+yp(3)+hs(1))/2-hs(1)/2 ws(3) hs(1)]);
    hold off;
    h=plot(x(pid)-0.01,erfs(pid,1),'g');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');
    hold on;
    h=plot(x(pid)+0.01,erfs(pid,1),'r');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');

    nl = x.^2;
    subplot('Position',[xp(3)+ws(3)  dy+0.05+(yp(2)+yp(3)+hs(1))/2-hs(1)/2 ws(2) hs(2)]);
    h = plot(nl,'k');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');

    dy = 0.4;
    subplot('Position',[xp(3) dy+(yp(2)+yp(3)+hs(1))/2-hs(1)/2 ws(3) hs(1)]);
    hold off;
    h=plot(x(pid)-0.01,orfs(pid,1),'g');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');
    hold on;
    h=plot(x(pid)+0.01,orfs(pid,1),'r');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');

    nl = x.^2;
    subplot('Position',[xp(3)+ws(3)  dy+0.05+(yp(2)+yp(3)+hs(1))/2-hs(1)/2 ws(2) hs(2)]);
    h = plot(nl,'k');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');



    dy = -0.3;
    subplot('Position',[xp(3) dy+(yp(2)+yp(3)+hs(1))/2-hs(1)/2 ws(3) hs(1)]);
    hold off;
    h=plot(x(pid)-0.01,erfs(pid+20,1),'g');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');
    hold on;
    h=plot(x(pid)+0.01,erfs(pid+20,1),'r');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');

    dy = -0.3;
    dx=0.25;
    subplot('Position',[dx+xp(3) dy+(yp(2)+yp(3)+hs(1))/2-hs(1)/2 ws(3) hs(1)]);
    hold off;
    h=plot(x(pid)-0.01,erfs(pid,2),'g');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');
    hold on;
    h=plot(x(pid)+0.01,-erfs(pid,2),'r');
    set(h,'linewidth',2); axis('tight'); set(gca,'xtick',[],'ytick',[],'color','none');



end
