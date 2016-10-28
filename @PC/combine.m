function AllExpt = combine(DATA, varargin)
%Try calling combine functions from PlotClusters
strargs = cell2cellstr(varargin);

Expts = getappdata(DATA.toplevel,'Expts');
D.Expts = Expts;
D.expid = find(strcmp(DATA.markexpts, DATA.expnames));
D.combineids = D.expid;
D.state.online = 0;
D.extype = 2;

D.state.useguilist = 0;
cmb.CombinePlot(D,0,'cluster',1,'ids',D.expid);

if sum(strcmp('allcells',strargs))
    C = PC.GetValue(DATA, 'CellList');
    cells = unique(C(D.expid,:,:));
    cells = cells(cells>0);
    for j = 1:length(cells)
        a = find(C == cells(j));
        [e,p,cl] = ind2sub(size(C),a);
        id = find(ismember(e,D.expid));
    end
    cmb.CombinePlot(D,0,'cluster',1,'ids',e(id));
end