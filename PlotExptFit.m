function h = PlotExptFit(fit, varargin)
%PlotExptFit(fit) Plots result from FitExpt
%if fit is a cell array, superimposes fits
%...,'label') places label at peak of fit

color = [];
matchstr = [];
showlabel = '';
colorset = 0;

if iscell(fit)
    for j = 1:length(fit)
        h{j} = PlotExptFit(fit{j}, varargin{:});
    end
    return;
end

if length(fit) > 1
    colors = mycolors;
    for j = 1:length(fit)
        h{j} = PlotExptFit(fit(j), 'color', colors{j},varargin{:});
    end
    return;
end

j = 1;
for j = 1:length(varargin)
    if strncmpi(varargin{j},'color',5)
        j = j+1;
        color = varargin{j};
        colorset = 1;
    elseif strncmpi(varargin{j},'match',5)
        j = j+1;
        matchstr = varargin{j};
    elseif strncmpi(varargin{j},'label',5)
        showlabel = 'name';
    end
    j = j+1;
end

if ~isempty(matchstr)
    if isfield(fit,'expt') & regexp(fit.expt,'matchstr')
        go = 1;
    else
        return;
    end
end

if isfield(fit, 'xfit')
    if isempty(color)
    color = mycolors;
    end
    for k = 1:size(fit.resp,2)
        h(k) = plot(fit.xfit,fit.fitted(:,k),'color',color{k});
        labels{k} = sprintf('%s=%f',fit.type{2},fit.y(k));        
        hold on;
    end
    mylegend(h,labels);
else
    h(1) = plot(fit.xv,fit.fitcurve);
end
if isfield(fit,'resp')
    hold on;
    sem = [];
    if isfield(fit,'sem')
        sem = fit.sem;
        for k = 1:size(fit.resp,2)
            h(end+1) = errorbar(fit.x, fit.resp(:,k),sem(:,k),'o');
            set(h(end),'color',color{k},'markerfacecolor',color{k});
        end
    elseif isfield(fit.state,'sds')
        if isfield(fit.state,'nreps')
            sem = fit.state.sds./sqrt(fit.state.nreps);
        else
            sem = fit.state.sds;
        end
        h(2) = errorbar(fit.x, fit.resp,sem,'o');
    end
    if isempty(sem)
        a = plot(fit.x, fit.resp,'o');
        h(end+1:end+length(a)) = a;
    end
end
if colorset
    set(h,'color',color);
end
if ~isempty(showlabel) && isfield(fit,'name')
    [a,b] = max(fit.fitcurve);
    [c,filename] = fileparts(fit.name);
    if isfield(fit,'probe')
        s = sprintf('%sP%d',filename,fit.probe);
    else
        s = sprintf('%s',filename);
    end
    text(fit.xv(b),a,s);
end

if isfield(fit,'extras')
    ex = fit.extras;
    sem = ex.sd./sqrt(ex.n);
   for j = 1:length(ex.label)
       h(end+1) = errorbar(ex.x(j),ex.means(j),sem(j),'o');
       text(ex.x(j),ex.means(j),ex.label{j},'horizontalalignment','left');
   end
end

if isfield(fit,'loadname')
    title(GetName(fit.loadname,'withcell','withexpt'));
elseif isfield(fit,'name')
    if isfield(fit,'probe')
        title(sprintf('%sP%d',GetName(fit.name),fit.probe));
    elseif isfield(fit,'cell')
        title(sprintf('%s Cell%d',GetName(fit.name),fit.cell));
    else
        title(fit.name);
    end
end