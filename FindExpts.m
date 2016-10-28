function list = FindExpts(path, type, varargin)
%list = FindExpts(path, type, ...) find data files that contain expts
%1matching type. path can be a cell array of strings
%
%FindExpts(list, type) will search within previous results.  e.g.  
%
%  list = FindExpts({'/b/data/lem/M*' '/b/data/jbe/M*'},'rls');
%  FindExpts(list, 'rls.*FrRC');
%
%step 1 builds list of all expts with 'rls' in the expt name
%step 2 extract just those that match the new pattern
%
%  FindExpts(list, pattern); prints total number of expts run that match 
%                   (NB can have > 1 cell per expt)
%  FindExpts(list, pattern,'list'); prints a list of sessions and numbers of trials
%  FindExpts(list, pattern,'listfile'); Checks each location for combined
%                data files, and prints out a list, with numbers of trials and cells
%          This requires loading the files, so can take a few minutes. It
%          returns a structure with summary data that can then be probed
%          further:
%
%  x = FindExpts(list, pattern,'listfile') 
%  FindExpts(x,'list') prints list of locations, with # cells found
%  FindExpts(x,'cells') same list, but locaionts with cells listed first
%  FindExpts(x,'cellsonly') just list locations with combined data found
%  FindExpts(x,'missing') prints list of expts run that have no combined files
%                    also prints out how many trials were saved
%
%   FindExpts(x,'celllist') prints trial count and spikerate for each cell
%
%  data = FindExpts(list, pattern,'listfile','fit') makes summary struct
%runs FitExpt on all cell data found, and records simple measures -
%respvar for   subspace data
%DDI/anova for tuning curves
%P = FindExpts(data,'celllist') generates a plot with tuning strength metrics
%        And returns a structure that can be used to refine the plot with
%        expt.PlotList(P)
%
%see also ListExpts, expt.PlotList


if iscell(path) && (nargin == 1 || sum(strcmp(type,{'cells' 'list' 'cellsonly' 'missing' 'celllist'})))
    if nargin == 1
        type =  'list';
    end
    DATA = [];
    for j = 1:length(path)
        DATA = PrintSummary(path{j},type,DATA,varargin{:});
    end
    if strcmp(type,'celllist') && DATA.n > 0 %extract DATA
        expt.PlotList(DATA, varargin{:});
    end
    list = DATA;
    return;
elseif iscell(path)
    for j = 1:length(path)
        [names{j}, details{j}] = ListExpts(path{j},type, 'silent', varargin{:});
    end
    allnames = cat(2, names{:});
    alldetails = cat(2, details(:));% ? redundant.  was details{:} but this makes struct, which is bad
    list = ListExpts(alldetails, type); % find expt
%    list = ListExpts(list,type); %find cells in file % not needed any more
    return;
elseif isfield(path,'allnames') && isfield(path,'cells')
    for j = 1:length(path.matchnames)
        for k = 1:length(path.matchnames{j})
            if path.ncells(j) == 0
            fprintf('%s:%s %d Trials No cells\n',path.filename{j}{k},path.matchnames{j}{k},path.ntrials(j));
            else
                for c = 1:length(path.cells{j})
                    fprintf('%s:%s Cell %d %d Trials\n',...
                        path.filename{j}{k},path.matchnames{j}{k},...
                        path.cells{j}(c),path.celltrials{j}(c));
                end
            end
        end
    end
    return;
elseif isfield(path,'allnames')
    list = ListExptTypes(path,type,varargin{:});
    return;
end


[names, details] = ListExpts(path,'silent');
if ~isempty(details)
    list = ListExpts(details, type);
    ListExptTypes(list);
end



function datalist = ListExptTypes(list, type, varargin)
j = 1;
strargs = cell2cellstr(varargin);
recalc = cellstrcmp('recalc',varargin);
    matchnames = {};
    for j = 1:length(list.matchnames)
        if nargin > 1
            id = find(CellToMat(regexp(list.matchnames{j},type)));
            matches{j} = id;
            matchnames = cat(2,matchnames, list.matchnames{j}(id));
        else
            matchnames = cat(2,matchnames, list.matchnames{j});
            matches{j} = 1:length(list.matchnames{j});
            matchlist{j} = ones(size(matches{j})) .* j;
        end
    end
    [exn, extypes] = Counts(matchnames);
    if isempty(exn)
        fprintf('No Expts match %s\n',type);
        datalist = {};
    elseif sum(strncmp('listfile',strargs,8))
        for j = 1:length(extypes)
            datalist{j} = ListExpts(list,extypes{j},'exact', varargin{:});
        end
        DATA.n = 0;
        for j = 1:length(datalist)
            DATA = PrintSummary(datalist{j},DATA);
        end
    elseif sum(strncmp('list',strargs,4))
        n = 0;
        for j = 1:length(exn)
%            fprintf('%s %d expts\n',extypes{j},exn(j));
            datalist{j}.expname = extypes{j};
            for k = 1:length(list.matchnames)
                id = find(strcmp(extypes{j},list.matchnames{k}));
                for e = id(:)'
                    if e > length(list.filename{k})
                        name = list.filename{k}{end};
                    else
                        name = list.filename{k}{e};
                    end
                    if strcmp(extypes{j},'rds.dxXce') & list.n2(k,e) > 3
                        fprintf('*');
                    end
                    fprintf('%s %s %d trials run %dX%d stim\n',name,list.matchnames{k}{e},list.trialcounts(k,e),list.nt(k,e),list.n2(k,e));
                end
                if ~isempty(id)
                    n = n+1;
                    datalist{j}.names{n} = [name ' ' list.matchnames{k}{id(1)}];
                    datalist{j}.nblocks(n) = list.matchcounts(k,id(1));
                    datalist{j}.ntrials(n) = list.trialcounts(k,id(1));
                    datalist{j}.dirpath{n} = list.dirpath{k};
                    datalist{j}.filename{n} = list.matchnames{k}{e};
                    if length(id) > 1 %for debugger, should not happen
                        fprintf('TWo Expts Match %s in %s\n',extypes{j},list.filename{k}{e});
                        datalist{j}.nblocks(n) = sum(list.matchcounts(k,id(1)));
                    end
                end                
            end
        end
    else
        allmatchlist = cat(2,matchlist{:});
        for j = 1:length(exn)
            datalist{j}.expname = extypes{j};
            datalist{j}.nexpts = exn(j);
            id = find(strcmp(extypes{j},matchnames));
            id = allmatchlist(id);
            if exn(j) > 1
                fprintf('%s %d expts\n',extypes{j},exn(j));
            else
                fprintf('%s %d expt in %s %d trials\n',extypes{j},exn(j),list.dirpath{id},sum(list.trialcounts(id,:)));
            end
        end
    end
    
    function DATA = PrintSummary(D, varargin)
        strargs = cell2cellstr(varargin);
recalc = cellstrcmp('recalc',varargin);
        DATA.n = 0;
        searchname = {};
        j = 1;
        while j <= length(varargin)
            if isfield(varargin{j},'n')
                DATA = varargin{j};
            elseif strcmp(varargin{j},'name')
                j = 1+j;
                searchname = {searchname{:} varargin{j}};
            end
            j = j+1;
        end
        
        if ~isempty(searchname)
            good = zeros(size(D.dirpath));
            id = find(CellToMat(regexp(D.dirpath,searchname{1})));
            good(id) = 1;
        else
            good = ones(size(D.dirpath));
        end
        if sum(strcmp('celllist',strargs))
            id = find(D.ncells > 0 & good);
            DATA.crit.nt = 1;
            if ~isfield(D,'respvar') && ~isfield(D,'di')
                fprintf('Need to Run FindExpts with ''listfile'',''fit'' before Using this option\n');
                return;
            end
            for j = id(:)'
                for k = D.cells{j}(:)'
                    ok = 0;
                    if isfield(D,'respvar') && length(D.respvar{j}) >= k
                        xstr = sprintf(' RespVar %.2f %.2f',D.respvar{j}(k).respvar,D.respvar{j}(k).skew);
                        if isempty(D.respvar{j}(k).respvar)
                            fprintf('No data for %s Cell %d\n',D.files{j}{k},D.cells{j}(k));
                        else
                            DATA.n = DATA.n+1;
                            ok = 1;
                            if recalc
                                [x, xx] = rc.RespVar(D.plots{j}{k});
                                D.respvar{j} = CatStruct(D.respvar{j},xx,k);
                            end
                            nv = length(D.respvar{j}(k).respvar);
                            DATA.respvar(DATA.n,1:nv) = D.respvar{j}(k).respvar;
                            DATA.Frs(DATA.n,1:nv) = D.respvar{j}(k).Frs;
                            DATA.ntrials(DATA.n,1:nv) = D.respvar{j}(k).n;
                            DATA.skew(DATA.n) = D.respvar{j}(k).skew;
                            if isfield(D.respvar{j},'respvarx')
                                DATA.respvarx(DATA.n,1:nv) = D.respvar{j}(k).respvarx;
                                ns = length(D.plots{j}{k}.sdfs.s(:));
                                DATA.nsdf(DATA.n) = ns;
                                DATA.pval(DATA.n,1:nv) = 1-fcdf(DATA.respvarx(DATA.n,1:nv),ns-1,ns-1);
                            end
                        end
                            
                    elseif isfield(D,'di')  && length(D.di{j}) >= k
                        DATA.n = DATA.n+1;
                        ok = 1;
                        DI = D.di{j}(k);
                        if isfield(DI,'di')
                            xstr = sprintf(' DI %.2f %.2f',DI.di,mean(DI.var));
                            DATA.di(DATA.n) = DI.di;
                            DATA.var(DATA.n) = mean(D.di{j}(k).var);
                        elseif isfield(DI,'ddi')
                            xstr = sprintf(' DI %.2f %.2f',DI.ddi,DI.sqerr);
                            DATA.di(DATA.n) = DI.ddi;
                            DATA.var(DATA.n) = DI.sqerr;
                            DATA.nt(DATA.n,1) = prctile(DI.ntrials,50);
                            if isfield(DI,'rates')
                                clear f;
                                for t = 1:length(DI.rates)
                                    f(t) = mean(DI.rates{t});
                                end
                                DATA.maxrate(DATA.n) = max(f);
                            else
                                DATA.maxrate(DATA.n) = NaN;
                            end
                        else
                            xstr = sprintf(' DI %.2f %.2f',DI.alldi,mean(D.di{j}(k).var));
                            DATA.di(DATA.n) = DI.alldi;
                            DATA.var(DATA.n) = mean(D.di{j}(k).var);
                        end
                        DATA.anovap(DATA.n) = D.di{j}(k).anovap;
                    else
                        xstr = '';
                    end
                    if isfield(D,'fits') && ok
                        DATA.fits{DATA.n} = D.fits{j}{k};
                    end
                    if isfield(D,'plots')
                        DATA.plots{DATA.n} = D.plots{j}{k};
                        R = D.plots{j}{k};
                        if isfield(R,'n')
                            DATA.nt(DATA.n,2) = prctile(R.n(R.n>0),50);
                            DATA.nt(DATA.n,3) = prctile(R.n(R.n>0),75);
                        end
                        if isfield(R,'Header')
                            x = mean(R.Header.dips,1);
                            if x(5) > 1000
                                DATA.isolation(DATA.n) = x(4);
                            elseif 1
                                DATA.isolation(DATA.n) = -x(1);
                            else
                                DATA.isolation(DATA.n) = max(x(2:3));
                            end
                            DATA.dropi(DATA.n) = mean(R.Header.dropi);
                        end
                    end
                    DATA.name{DATA.n} = D.cellfiles{j}{k};
                    DATA.monkeyname{DATA.n} = GetMonkeyName(D.cellfiles{j}{k});
                    DATA.cellnumber(DATA.n) = k;
                    fprintf('%s/%s Cell %d/%d %d trials at %.1f spk/sec%s\n',D.dirpath{j},D.filename{j},k,D.ncells(j),D.celltrials{j}(k),D.rates{j}(k),xstr);
                end
            end
        elseif sum(strncmp('cells',strargs,5))
            id = find(D.ncells > 0);
            for j = id(:)'
                fprintf('%s/%s %d cells %d files %d trials\n',D.dirpath{j},D.filename{j},D.ncells(j),D.nfiles(j),D.ntrials(j));
            end
            if sum(strcmp('cellsonly',strargs)) == 0
                id = find(D.ncells == 0);
                for j = id(:)'
                    fprintf('%s/%s %d cells %d files %d trials\n',D.dirpath{j},D.filename{j},D.ncells(j),D.nfiles(j),D.ntrials(j));
                end
            end
        elseif sum(strncmp('missing',strargs,5))
            if ~isfield(D,'ncells')
                fprintf('Cannot list missing cells unles you have run FindExpts(list,''listdata'')\n');
                return;
            end
                id = find(D.ncells == 0 & D.nblock > 0 & good);
                for j = id(:)'
                    fprintf('%s/%s %s %d cells in %d disk files %d trials. %d Blocks,%d Trials Run\n',D.dirpath{j},D.filename{j},D.pattern,D.ncells(j),D.nfiles(j),D.ntrials(j),D.nblock(j),D.trialsrun(j));
                    if isfield(D,'Comments') && length(D.Comments) >= j 
                        for k = 1:length(D.Comments{j})
                            cprintf('blue','%s:%s\n',D.dirpath{j},D.Comments{j}{k});
                        end
                    end
                end
        else
            for j = 1:length(D.dirpath)
                if isfield(D,'ncells') %did do search through completed files
                    fprintf('%s/%s %d cells %d files %d trials. %d Blocks run\n',D.dirpath{j},D.filename{j},D.ncells(j),D.nfiles(j),D.ntrials(j),D.nblock(j));
                else %only have list of what expts were dones
                    fprintf('%s/%s %d blocks run %d trials\n',D.dirpath{j},D.filename{j},D.nblocks(j),D.ntrials(j));
                end
            end
        end
        

            
