function PlotPlaid(res, varargin)
%PlotPlaid(res, varargin) plots subsets of plaid expt conditions
%default is to to show 1D, UniKinetic, and Type I in a polar plot

figtag = 'Plaid';
plottype = 'component';
vskip = [];
j = 1;
while j <=length(varargin)
    if strcmp(varargin{j},'tag')
        j = j+1;
        figtag = varargin{j};
    elseif strncmp(varargin{j},'withflicker',6)
        vskip = [vskip j];
        plottype = 'withflicker';
    elseif strncmp(varargin{j},'withrds',6)
        vskip = [vskip j];
        plottype = 'withrds';
    end
    j = j+1;
end
varargin = varargin(setdiff(1:length(varargin),vskip));

if iscell(res)
    for j = 1:length(res)
        name = GetName(res{j},'withcell');
        PlotPlaid(res{j},varargin{:},'tag',name);
    end
    return;
end
if isexpt(res)
    Expt = res;
    res = PlotExpt(Expt,'noplot');
end
  

C = res.conditions;
a2vals = GetCondition(C,'a2');
c2vals = GetCondition(C,'c2');
sM = GetCondition(C,'sM');
v2 = GetCondition(C,'v2');
sl = GetCondition(C,'sl');
st = GetCondition(C,'st');
oneid = find(c2vals == 0 & (sM == 0 & ~isinf(v2)) & sl ==1 & st ~= 2);
colors{4} = [0 0 0];
colors{2} = [1 0 0];
colors{3} = [0 0 1];
colors{1} = [0 1 0];
id = find(a2vals ~= 90 & c2vals == 1 & (sM == 0 & ~isinf(v2)));
flid = find(a2vals ~= 90 & c2vals == 1 & (sM == 33 | isinf(v2)));
rdsid = find(st ==2);
for j = 1:length(colors)
    linestyles{j} = '-';
end
facecolors = colors;

if strcmp(plottype,'withflicker')
    colors(5) = colors(1);
    colors(6) = colors(2);
    linestyles{5} = '--';
    linestyles{6} = '--';
    facecolors{5} = 'none';
    facecolors{6} = 'none';
    PlotResult(res,'polar',varargin{:},'bvals',[id' oneid' flid'],'tag',figtag,'colors',colors,facecolors,'linestyles',linestyles);
elseif strcmp(plottype,'withrds')
    colors(5) = colors(4);
    facecolors{5} = 'none';
    linestyles{5} = '--';
    PlotResult(res,'polar',varargin{:},'bvals',[id' oneid' rdsid'],'tag',figtag,'colors',colors,facecolors,'linestyles',linestyles);
else
    PlotResult(res,'polar',varargin{:},'bvals',[id' oneid'],'tag',figtag,'colors',colors);
end

