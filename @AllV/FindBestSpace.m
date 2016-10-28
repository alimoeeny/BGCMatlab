function [bestspace, DATA] = FindBestSpace(DATA, C, varargin)
%AllV.FindBestSpace(DATA, C, varargin)
%Given a Cluster classification, find an ND space that optimizes the
%separation measure using mahal dprime
% by default ccompares DATA.currentcluster vs 0
%AllV.FindBestSpace(DATA,C,'clusters',a,b)
% Finds space where mahal distance of b, from distribution defined by a, is
% max. I.e want cell first, hash second. 

if isfigure(DATA)
    DATA = get(DATA,'UserData');
    C = DATA.cluster;
end
   
cl{1} = DATA.currentcluster+1;
cl{2} = 0;
plottype = 'none';
verbose = 'quiet';
isolationtest = 'sumahal';
searchmode = 'quick';
method = 2;
usethree = 1;
subspace = [];
plotspace = [];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'allclusters',5) %look at all combinations
        [bestspace, DATA] = FindAllSpace(DATA,C);
        if nargout < 2
            SetData(DATA);
        end
        return;
    elseif strncmpi(varargin{j},'clusters',5) %these are id # in clst. So Clsuter 1 == 2
        j = j+1;
        cl{1} = varargin{j};
        j = j+1;
        cl{2} = varargin{j};        
    elseif strncmpi(varargin{j},'exhausttest',9)
        method = 4;
    elseif strncmpi(varargin{j},'exhaustive',5)
        method = 2;
    elseif strncmpi(varargin{j},'nearmahal',5)
        isolationtest = varargin{j};
    elseif strncmpi(varargin{j},'plotspace',5)
        j = j+1;
        plotspace = varargin{j};
    elseif strncmpi(varargin{j},'plotallxy',5)
        plottype = 'plotallxy';
    elseif strncmpi(varargin{j},'plotxy',5)
        plottype = 'plotxy';
    elseif strncmpi(varargin{j},'quick',5)
        method = 3;
    elseif strncmpi(varargin{j},'twod',4)
        usethree = 0;
    elseif strncmpi(varargin{j},'subspace',4)
        j = j+1;
        subspace = varargin{j};
    elseif strncmpi(varargin{j},'verbose',5)
        verbose = varargin{j};
    end
    j = j+1;
end

idstr = [AllV.IDStr(DATA) ':'];
if DATA.usegmcid
    gid = find(ismember(DATA.gmcid,cl{1}));
    nid = find(ismember(DATA.gmcid,cl{2}));
else
    gid = find(ismember(DATA.clst,cl{1}));
    nid = find(ismember(DATA.clst,cl{2}));
    if cl{2} == 0
        nid = find(~ismember(DATA.clst,cl{1}));
    end        
end

if isfield(C,'isolation')
    isolation = C.isolation;
else
    isolation = [NaN NaN];
end

%this determine which cluster defines the distribution
%CalcIsoltion uses distance of 1 from distrbution given by 2
%cl{1} is gid, so cl{1} is treated as the SU
    pts = [nid(:)' gid(:)'];
    clst = DATA.clst;
    clst(nid) = 1;
    clst(gid) = 2;

if isfield(C,'space')
    if C.space(1) == AllV.USERSPACE %klugde for now. 7 is arbitray space, usually TemplateScores
        C.space(1) = AllV.TEMPLATESPACE;
    end
    bestspace.space = C.space;
end

T = [];
O = [];
oa = 0;
oscore = [0 0 0 0];
altspace = '';
realspace = C.space(1);
if isfield(DATA,'spikefeatures') && strcmpi(C.autocutmode,'james')
    T = DATA.spikefeatures;
    C.space = 10;
    for j = 1:size(T,2)
        Labels{j} = sprintf('Spike Feature(James) %d',j);
    end
elseif C.space(1) == 3
    T = DATA.TemplateScores;
    Labels = DATA.TemplateLabels;
    if method ==2 && isappdata(DATA.toplevel,'AllTemplates')
        A = getappdata(DATA.toplevel,'AllTemplates');
        if length(A) > 1
            T = cat(2,A(1).Scores,A(2).Scores);
            for j = 1:length(DATA.TemplateLabels)
                Labels{j} = sprintf('T1.%s',DATA.TemplateLabels{j});
                Labels{j+length(DATA.TemplateLabels)} = sprintf('T2.%s',DATA.TemplateLabels{j});
            end
            realspace = AllV.USERSPACE;
        end
    end
    if isempty(subspace)
        [~, details] = RemoveRedundantRows(T,0.9,idstr);
        subspace = details.good;
    else
        [~, details] = RemoveRedundantRows(T(:,subspace),0.9,idstr);
        subspace = subspace(details.good);
    end
    O = DATA.pcs;
    if size(O,2) > 10
        O = O(:,1:6);
    end
    altspace = 'PC';
    
elseif C.space(1) == 1
    T = DATA.pcs;
    if size(T,2) > 10
        T = T(:,1:6);
    end
    for j = 1:size(T,2)
        Labels{j} = sprintf('PC %d',j);
    end

elseif C.space == 2
    T = [];
else
    T = [];
end

bestspace.SU = cl{1};
bestspace.MU = cl{2};
bestspace.xnames = {};
if ~isempty(T)
    subspace = union(subspace,plotspace);%make sure dimensions for plotting are kept
    if ~isempty(subspace)
        T = T(:,subspace);
        plotspace = find(ismember(subspace,plotspace));
    else
        subspace = 1:size(T,2);
    end
    bestspace.ndims = size(T,2);
    if method ==1 %just look at this cluster vs everyone else
        dims = C.space(2:3);
        newdims = setdiff(1:size(T,2),dims);
        d(1,:) = CalcIsolation(T(pts,dims),clst(pts),2);
        for j = 1:length(newdims)
            d(j+1,:) = CalcIsolation(T(pts,[dims newdims(j)]),clst(pts),2);
        end
        max(d);
    elseif method == 3
        [x, dims, ddetails] = clust.FindBestSpace(T(pts,:),clst(pts),2, isolationtest);
        bestspace.isolation = ddetails.isolation;

    elseif method ==2 || method == 4 %try all pairings of spaces
        ts = now;
        if method == 4            
            [x, xdims, ddetails] = clust.FindBestSpace(T(pts,:),clst(pts),2, isolationtest);
            bestspace.quickisolation = ddetails.isolation;
            bestspace.quickspace = xdims;
            bestspace.dur(1) = mytoc(ts);
        end
        n = 1;
        for j = 2:size(T,2)
            for k = 1:j-1;
                [score(n,:), details] = CalcIsolation(T(pts,[j k]),clst(pts),2);
                space(n,1) = j;
                space(n,2) = k;
                n = n+1;
            end
        end
        for j = 2:size(O,2)
            for k = 1:j-1;
                oscore(n,:) = CalcIsolation(O(pts,[j k]),clst(pts),2);
                ospace(n,1) = j;
                ospace(n,2) = k;
                n = n+1;
            end
        end
        sid = [];
        xdims = [];
        if size(score,1) > 1
            [a,b] = max(score);
            [~,sid] = sort(score(:,1),'descend');
        else
            a = score;
            b =1;
        end
        if strcmp(isolationtest,'nearmahal')
            dims = space(b(3),:);
        elseif strcmp(isolationtest,'mahal') || strcmp(isolationtest,'sumahal')
            dims = space(b(1),:);
            if ~isempty(sid)
                xdims = space(sid(1:min([length(sid) 5])),:);
            end
        else
            dims = unique(space(b,:))';
            if length(dims) > 2
                d = CalcIsolation(T(pts,dims),clst(pts),2);
                if sum(d-a) < 0.5  %no real improvement
                    dims = space(b(1),:);
                else
                    a = d;
                end
            end
        end
        newdims = setdiff(1:size(T,2),dims);
        ds(1,:) = CalcIsolation(T(pts,dims),clst(pts),2);
%
        if usethree && size(T,2) > 2
            for j = 1:length(newdims)
                ds(j+1,:) = CalcIsolation(T(pts,[dims newdims(j)]),clst(pts),2);
            end
        [c,d] = max(ds); %max of mutual distances, not hte "near" calc
        if strcmp(isolationtest,'nearmahal')
            if c(3) - a(3) > 0.1 && d(3) > 1
                dims = [dims newdims(d(3)-1)];
                bestspace.isolation = ds(d(3),:);
                bestspace.twodisolation = ds(1,:);
            else
                bestspace.isolation = ds(1,:);
            end
        else
            if sum(c(1:2)-a(1:2)) > 0.5 %3D space really better
                if c(1) - a(1) > 0.1 %in case c() ~ a(1) and c(2) >> a(2)
                    if d(1) > 1
                        dims = [dims newdims(d(1)-1)];
                    else
                        cprintf('red','BestSpace: 3D > 2D but no third dimension used!\n');
                    end
                    bestspace.isolation = ds(d(1),:);
                else
                    if d(2) > 1
                    dims = [dims newdims(d(2)-1)];
                    bestspace.isolation = ds(d(2),:);
                    bestspace.twodisolation = ds(1,:);
                    else
                        cprintf('red','!!Error - dimensions wrong for best space\n');
                    end
                end
            else
                bestspace.isolation = ds(1,:);
                if sum(d > 1)  %a 3d space was a bit better
                    dm =[ 0 0 ];
                    id = find(d > 1);
                    dm(id) = newdims(d(id)-1);
                    if strcmp(verbose,'verbose')
                        fprintf('Best 3D space isolation %.3f,%.3f including space %d,%d\n',ds(d(1),1),ds(d(2),2),dm(1),dm(2));
                    end
                end
            end
        end
        else
            bestspace.isolation = ds;
        end
        if isfield(bestspace,'dur')
            bestspace.dur(2) = mytoc(ts) - bestspace.dur(1);
        end
        [oa,ob] = max(oscore);
    end
    if ~isempty(xdims)
        xdims = setdiff(unique(xdims),dims);
        xdims = subspace(xdims);
    end
%plot before reversing subspace    
    if ~isempty(plotspace)
        AllV.SetFigure(DATA.tag.tmplscore,DATA);
        if 0
            PlotND(T(:,plotspace),[],'spkid',DATA.clst);
            xlabel(Labels{plotspace(1)});
            ylabel(Labels{plotspace(2)});
        else
            [Ev,E] = eig(cov(T(:,plotspace)));
            E = diag(E);
            [~,order] = sort(E,'descend');
            pcs = T(:,plotspace) *Ev(:,order(1:2));
            PlotND(pcs,[],'spkid',DATA.clst);
            DATA.TemplatePC.dims = subspace(plotspace);
            DATA.TemplatePC.labels = Labels(subspace(plotspace));
            DATA.TemplatePC.ev = Ev(:,order);
            eva = sum(abs(Ev(:,1:2)),2);
            [~,id] = sort(eva,'descend');
            DATA.TemplatePC.dimid = id;
            for j = 1:length(id)
                fprintf('Space %s weight %.2f\n',DATA.TemplatePC.labels{id(j)},eva(id(j)));
            end
            setappdata(gcf,'ClusterSpace',[8 1 2]);
            X.Variables{1} = 'TemplatePC1';
            X.Variables{2} = 'TemplatePC2';
            xlabel(X.Variables{1});
            ylabel(X.Variables{2});
            X.pcplot = [1 2];
            X.ClusterSpace = [8 12];
            set(gca,'UserData',X);
            if nargout < 2
                SetData(DATA);
            end
        end
    else
        if strcmp(plottype,'plotxy')
            AllV.SetFigure(DATA.tag.tmplscore,DATA);
            PlotND(T(pts,dims),[],'spkid',clst(pts));
        end
        if strcmp(plottype,'plotallxy')
            AllV.SetFigure(DATA.tag.tmplscore,DATA);
            PlotND(T(:,dims),[],'spkid',DATA.clst);
        end
        xlabel(Labels{subspace(dims(1))});
        ylabel(Labels{subspace(dims(2))});
    end
    if ~isempty(subspace)
        dims = subspace(dims);
    end
    bestspace.space = [C.space(1) dims(:)'];
    if C.space(1) == 3
        xstr = sprintf('(%s vs %s)',Labels{dims(1)},Labels{dims(2)});
    else
        xstr = '';
    end
    bestspace.str =xstr;
    if strcmp(verbose,'verbose')
        fprintf('Best Mahal distance cl%d vs cl%d  %.3f,%.3f in space%s%s (%.3f,%.3f native)\n',cl{1}-1,cl{2}-1,bestspace.isolation,sprintf(' %d', bestspace.space),xstr,isolation(1),isolation(2));
        if oa(1) > bestspace.isolation(1)
            fprintf('Better score in %s space: %.2f',altspace,oa(1));
        end
    end
    if realspace == AllV.USERSPACE
        bestspace.space(1) = realspace;
        bestspace.Variables{1} = Labels{dims(1)};
        bestspace.Variables{2} = Labels{dims(2)};
    end
    if ~isempty(xdims)
        bestspace.xdims = xdims;
        bestspace.xnames = Labels(xdims);
    end
end


function [space, DATA] = FindAllSpace(DATA, C, varargin)

cls = unique(DATA.clst);
allspace = [];
spaces = {};
for j = 2:length(cls)
    for k = 1:j-1;
        fprintf('Checking %d vs %d\n',j,k);
        spaces{end+1} = AllV.FindBestSpace(DATA,C,'clusters',cls(j),cls(k));
        allspace = [allspace spaces{end}.space(2:end)];
        if k > 1 %both are cells, try both ways
            spaces{end+1} = AllV.FindBestSpace(DATA,C,'clusters',cls(k),cls(j));
            allspace = [allspace spaces{end}.space(2:end)];
        end
    end
end
allspace = unique(allspace);
[space, DATA] = AllV.FindBestSpace(DATA,C, 'plotspace', allspace);
space.space = allspace;


