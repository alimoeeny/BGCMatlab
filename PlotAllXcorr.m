function [cells, xcorrs] = PlotAllXCorr(DATA, DataClusters, cells, varargin)
F = gcf;

callback = @PlotXcorr;
j = 1;
ids = [];
plotxc =1;
checktime = 1;
xctimes = -0.2005:0.001:0.2005;
efficonly = 0;
if isfield(DATA.plot,'showxc') && DATA.plot.showxc == 1
    plotxc = 0;
end


while j <= length(varargin)
    if strcmp(varargin{j},'callback')
        j = j+1;
        callback = varargin{j};
    elseif strcmp(varargin{j},'efficonly')
        xctimes = [-0.001:0.001:0.001];
        efficonly = 1;
        plotxc = 0;
    elseif strcmp(varargin{j},'sublist')
        j = j+1;
        ids = varargin{j};
    end
    j = j+1;
end
ts = now;

if isempty(ids)
    ids = 1:length(cells);
end
xcorrs = [];
    ClearPlot;
    if isfield(DATA,'toplevel') && isfigure(DATA.toplevel)
        setappdata(DATA.toplevel,'xcCellList',cells);
    end
%    xpts = linspace(DATA.spts(1),DATA.spts(end),401);
    nc = length(cells);
    nxc = 0;
    nspk = 0;
    exptno = 0;
    for j = 1:length(cells)
        cl = cells(j).cl -1;
        P = DataClusters{cells(j).p};
        if isfield(P,'exptno')
            exptno = P.exptno;
        end
        if cl > 0
            if length(P.next) >= cl && ~isempty(P.next{cl})
                P  = P.next{cells(j).cl-1};
                if ~isfield(P,'times')
                    P.times = [];
                end
            else
                P.times = [];
            end
        elseif isfield(P,'clst')
            if isfield(P,'t')
                DataClusters{cells(j).p}.times = P.t(P.clst==2);
            elseif isfield(P,'times') 
                if length(P.clst) <= length(P.times)
                    DataClusters{cells(j).p}.times = P.times(P.clst==2);
                else
                    cprintf('red','Expt %d Building Probe %d cluster %dfrom clst. Only have short  Times\n',exptno,cells(j).p,cells(j).cl);
                end
            else
                cprintf('red','Expt %d Building Probe %d cluster %d Missing Times\n',exptno,cells(j).p,cells(j).cl);
                DataClusters{cells(j).p}.times = [];
            end
        end
        spkcount(j) = length(P.times);
        nspk = nspk + length(P.times);
    end
    if checktime
        fprintf('Expt %d Building Xcorrs for %d (max(%d)) Total Spikes in %d cells\n',exptno,nspk,max(spkcount),nc);
    end
    if max(spkcount) > 10000
    end
    nspk = 0;
    allnspk = [0;0];
    lastts = now;
    for j = 1:length(cells)
        cl = cells(j).cl -1;
        P = PC.GetClusterInfo( DataClusters{cells(j).p},cells(j).cl);
        if checktime && mytoc(lastts) > 5
            fprintf('%d cells took%.2f\n',j,mytoc(ts));
            lastts = now;
        end
    for k = 1:j
        if ismember(k,ids) || ismember(j,ids)            
%        xoff = floor((j-1)/6) .* length(DATA.spts);
%        yoff = rem(j-1,6) .* DATA.vsep;

        Q = PC.GetClusterInfo( DataClusters{cells(k).p},cells(k).cl);
        nspk = [length(P.times) length(Q.times)];
        if prod(nspk) > 1000000 %watch out for memory issues
        end
        nxc = nxc+1;
        if efficonly
            [details.efficacy, edetails] = CalcEfficacy(P.times,Q.times);
            xcorrs(nxc).delays(1) = edetails.delay;
            if isfield(edetails,'xc')
               xc = edetails.xc;
               details.midpt = edetails.midpt;
            else
                xc = NaN;
            end
        else
            [xc, details] = xcorrtimes(P.times,Q.times,'times',xctimes);
            xpts = details.xpts;
        end
        xcorrs(nxc).nspk = nspk;
        xcorrs(nxc).shapexc = PC.ShapeCorr(P,Q,'delays',5);
        xcorrs(nxc).efficacy = details.efficacy;
        if isfield(DATA,'ArrayConfig')
            xcorrs(nxc).separation = ArrayDistance(DATA.ArrayConfig, cells(j).p,cells(k).p);
        else
            xcorrs(nxc).separation = NaN; %?use -abs(diff(p))???
        end

        if (j == k) && ~isnan(xc(1));
            xc(details.midpt) = 0;
        end
        xcorrs(nxc).xc = xc;
        xcorrs(nxc).p = [j k];
        xcorrs(nxc).probe = [cells(j).p + cells(j).cl/10  cells(k).p + cells(k).cl/10 ];
        
        if plotxc
            set(0,'CurrentFigure',F);
            mysubplot(nc,nc,(j-1) * nc+k);
            h = plot(xpts,xc,'k','linewidth',2);
            [a,b]= max(xc);
            if max(details.efficacy) > DATA.crit.synci(2)
                set(h,'color','r');
            end
            set(gca,'xtick',[],'ytick',[],'buttondownfcn',{callback, j,k});
            set(gca,'UserData',[j k]);
            set(h,'buttondownfcn',{callback, j,k});
            axis('tight');
            drawnow;
            if k == j
                if isfield(P,'fitdprime') && P.fitdprime(1) < -2
                    set(h,'color','r');
                else
                    set(h,'color','k');
                end
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                h = text(xl(2),yl(2),sprintf('%d/%d',cells(j).p,cells(j).cl),...
                    'horizontalalignment','right',...
                    'verticalalignment','bottom');
                if j == 1
                    set(h,'horizontalalignment','right',...
                        'verticalalignment','bottom');
                end
            elseif k == max(ids)
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                h = text(xl(2),yl(2),sprintf('%d/%d',cells(j).p,cells(j).cl),...
                    'horizontalalignment','left',...
                    'verticalalignment','top');
            elseif j == max(ids)
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                h = text(xl(2),yl(2),sprintf('%d/%d',cells(k).p,cells(k).cl),...
                    'horizontalalignment','right',...
                    'verticalalignment','bottom');
            end
        end
        xcorrs(nxc).ts = mytoc(ts);
        end
    end
    end
    mytoc(ts);
    if isfield(DATA,'toplevel') && isfigure(DATA.toplevel)
        setappdata(DATA.toplevel,'xcorrs',xcorrs);
    end
%    ReplotXcorrs(DATA, [], 'Shape/Efficacy')
    
function ReplotXcorrs(a,b, type)
    DATA = GetDataFromFig(a);
    xcorrs = getappdata(DATA.toplevel,'xcorrs');
    if sum(strcmpi(type, {'meanim' 'xcorrs' 'syncspikes' 'histograms'}))
        DATA.plot.xcorrtype = type;
        SetMenuCheck(a,'exclusive');
    elseif strcmpi(type, 'Shape/Efficacy')
        GetFigure(DATA.tag.xcorrpop);
        ClearPlot;
        hold off;
    for j = 1:length(xcorrs)
    plot(xcorrs(j).shapexc, max(xcorrs(j).efficacy),'o',...
        'buttondownfcn',{@PlotXcorr, xcorrs(j).p(1),xcorrs(j).p(2)});
    hold on;
    end
    
    set(gca,'yscale','log','ylim',[min([xcorrs.shapexc]) 1]);
    end

    
    function xc = xShapeCorr(P,Q)
    xc = corrcoef(P.MeanSpike.ms(:),Q.MeanSpike.ms(:));
    xc = xc(1,2);
