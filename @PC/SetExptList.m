function [DATA, Expts] = SetExptList(DATA, AllExpts, varargin)
%Expts = cmbSetExptList(DATA, AllExpts) Clean up Expts List, and align with
%Clusters
%DATA.exptids
%DATA.exptlist is list matching DATA.CellDetails.exptids (or equivalent)

for j = 1:length(AllExpts)
    if ~isempty(AllExpts{j})
        Expts{AllExpts{j}.Header.exptno} = AllExpts{j};
        goodexpt(j) = 1;
    end
end
%goodexpt = find(goodexpt ==1);
%     goodexpt = intersect(goodexpt,DATA.exptid); %Both Expt and ClusterTimes

%DATA.expid is list of expt #s for which ClusterTimes Data Exists.
%DATA.expnames{e} is name for expt DATA.exptid(e)
if length(goodexpt)
    for j = 1:length(Expts)
        if goodexpt(j)
            DATA.expnames{j} = Expt2Name(Expts{j},'addsuff');
            DATA.electrodedepth(j) = mean(GetEval(Expts{j},'ed'));
        else
            DATA.expnames{j} = '';
            DATA.electrodedepth(j) = NaN;
        end
    end
    setappdata(DATA.toplevel,'Expts',Expts);
    ExptList.expnames = DATA.expnames;
else
    ExptList.expnames = [];
end

it = findobj(DATA.toplevel,'tag','ExptList');
if ~isempty(it)
    str = get(it,'string');
    if length(Expts) > length(str)
        str = cellstr(num2str([1:length(Expts)]'));
        set(it,'string',str);
    end
end
ExptList.exptid = DATA.exptid;
if ~isfield(DATA,'exptlist')
    DATA.exptlist = DATA.exptid;
end
ExptList.usealltrials = DATA.usealltrials;
DATA.exptwithcluster = DATA.exptid;
if ~isempty(ExptList.expnames) %can happen with relist
    outname = [DATA.name '/ExptList.mat'];
    save(outname,'ExptList');
end
if nargout == 0
    SetData(DATA);
end
