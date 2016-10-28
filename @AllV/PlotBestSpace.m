function DATA = PlotBestSpace(DATA, varargin)
%Call PlotOneXY with space that shows bestspace

C = DATA.cluster;
cl = DATA.currentcluster;
isolation = [NaN NaN];
if ~isfield(C,'bestisolation')
    if C.auto && strcmp(C.autocutmode,'ecker')
        cl = DATA.currentcluster+1;
        if size(C.eckercluster.bestspace,1) < cl
            return;
        end
        space = C.eckercluster.bestspace(cl,:);
        isolation = C.isolation(1);
    end
else
    if cl > 1 && isfield(C.next{cl-1},'bestisolation')
        space = C.next{cl-1}.bestisolation.space;
        isolation = C.next{cl-1}.bestisolation.isolation;
        bestspace = C.next{cl-1}.bestisolation;
    else
        space = C.bestisolation.space;
        isolation = C.bestisolation.isolation;
        bestspace = C.bestisolation;
    end
end
AllV.SetFigure(DATA.tag.tmplscore, DATA);
if space(1) == 3 %template space
    xname = DATA.TemplateLabels{space(2)};
    yname = DATA.TemplateLabels{space(3)};
elseif space(1) == AllV.USERSPACE
    xname = bestspace.Variables{1};
    yname = bestspace.Variables{2};
    
elseif space(1) == 1
    xname = sprintf('PC%d',space(2));
    yname = sprintf('PC%d',space(3));
end
if exist('xname')
    DATA.xyplot.xy = {xname yname};
    xy = AllV.PlotOneXY(DATA,{ xname yname});
    SetData(DATA);
    title(sprintf('Isolation %.2f,%.2f',isolation));
    fprintf('Fitting Ellipse to classification\n');
    [E, score] = FindEllipse(xy, DATA.clst,'cluster',DATA.currentcluster);
    DrawEllipse(E,'add');
end