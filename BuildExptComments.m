function [Comments, ExptList, CellInfo] = BuildExptComments(path, varargin)
%BuildExptComments(path, varargin) Extracts comment lines from
% Spike2.mat files and puts them in path/ExptComments.mat
% By default does nothing in ExptComments.mat exists.
rebuild = 0;
showcomments = 0;
CellInfo = [];

j = 1;
while j <=length(varargin)
    if strncmpi(varargin{j},'rebuild',4)
        rebuild = 1;
    elseif strncmpi(varargin{j},'show',4)
        showcomments = 1;
    end
    j = j+1;
end

if ~isdir(path) && exist(path,'file')
    [path, filename] = fileparts(path);
end

outname = [path '/ExptComments.mat'];
if exist(outname) && ~rebuild
    load(outname);
else
    Comments = [];
    ExptList = [];
    if ~isdir(path)
        fprintf('No Folder %s\n',path);
        return;
    end

d = mydir([path '/*idx.mat']);
nc = 1;
ne = 1;
pes = [];
allxy = [];
for j = 1:length(d)
    clear Expt;
    load(d(j).name);
    if exist('Expt','var') && isfield(Expt,'Comments')
    Expt.Header.builddate = double(d(j).datenum);
    if isfield(Expt.Comments,'Peninfo') && isfield(Expt.Comments.Peninfo,'user')
        user = Expt.Comments.Peninfo.user;
    else
        user = '??';
    end
    for k = 1:length(Expt.Comments.text)
        if strncmp(Expt.Comments.text{k},' Electrode',7)
        elseif strncmp(Expt.Comments.text{k},'Electrode',7)
        elseif strncmp(Expt.Comments.text{k},'Experimenter Not Set',7)
        elseif strncmp(Expt.Comments.text{k},'rf=',3)
        elseif strncmp(Expt.Comments.text{k},'cm=NotSet',9)
        elseif strncmp(Expt.Comments.text{k},'cm=back=',8)
        else
            Comments(nc).text = Expt.Comments.text{k};
            Comments(nc).time = ConvertTime(Expt, Expt.Comments.times(k));
            Comments(nc).depth = FindElectrodeDepth(Expt, Expt.Comments.times(k));
            Comments(nc).user = user;
            nc = nc+1;
        end
    end
    if ~exist('ExptList','var') && isfield(Expt,'ExptList')
        Expts = ReadExptDir(path);
        for k = 1:length(Expts)
            ExptList(k).expname = Expt2Name(Expts{k});
            ExptList(k).start = Expts{k}.Header.trange(1);
            ExptList(k).end = Expts{k}.Header.trange(2);
        end
    end
    for k = 1:length(ExptList)
        ExptsList(ne).name = ExptList(k).expname;
        ExptsList(ne).start =  ConvertTime(Expt, ExptList(k).start);
        ExptsList(ne).end =  ConvertTime(Expt, ExptList(k).end);
        ExptsList(ne).depth(1) =  FindElectrodeDepth(Expt, ExptList(k).start);
        ExptsList(ne).depth(2) =  FindElectrodeDepth(Expt, ExptList(k).end);
        ExptsList(ne).suffix =  GetExptNumber(Expt);
        ExptsList(ne).ntrials =  CountTrials(Expt, [ExptList(k).start ExptList(k).end]);
        ne = ne+1;
    end
    if isfield(Expt.Trials,'rf') && size(Expt.Trials.rf,2) > 5
        pes(j) = median(Expt.Trials.rf(:,6));
        if size(Expt.Trials.rf,2) > 7
            allxy(j,1) = median(Expt.Trials.rf(:,7));
            allxy(j,2) = median(Expt.Trials.rf(:,8));
        end
    end
    end
end

if ~isempty(pes)
    pe = median(pes);
    penxy = median(allxy);
    logfile = sprintf('/b/bgc/anal/%s/pens/pen%d.log',GetMonkeyName(path),pe);
    if ~exist(logfile)
        logfile = sprintf('/b/data/%s/pens/pen%d.log',GetMonkeyName(path),pe);
    end
    
    pen = ReadPen(logfile,'noplot');
else
    pen = [];
    pe = 0;
    penxy = [NaN NaN];
end
if isempty(Comments)
    exptcomments = {};
else
    exptcomments = {Comments.text};
end
if isfield(pen, 'comments')
for j = 1:length(pen.comments)
    if ~strncmp(pen.comments{j},'NotSet',6) &&  isempty(strcmp(pen.comments{j},exptcomments))
        Comments(end+1).text = pen.comments{j};
        Comments(end).time = pen.date(pen.cmtime(j));
        Comments(end).depth = pen.depths(pen.cmtime(j));
    end
end
end


if ~isempty(Comments)
[a, b] = sort([Comments.time]);
Comments = Comments(b);
end
ExptList = ExptsList;
if ~isempty(ExptList)
[a, b] = sort([ExptList.start]);
ExptList = ExptList(b);
end
[a,b,c,d] = GetMonkeyName(path);
txtfile = [path '/' a c 'Add.txt'];
d = dir(txtfile);
if ~isempty(d)
txt = scanlines(txtfile);
id = find(strncmp('-2',txt,2));
for j = 1:length(id)
    Comments(end+1).text = txt{id(j)}(3:end);
    Comments(end).user = 'txt';
    Comments(end).time = d.datenum;
end
end
CellInfo.pen =[pe penxy];
save(outname,'Comments','ExptList','CellInfo');
end


function truet = ConvertTime(Expt, t)
%t in seconds

t = double(t);
if isfield(Expt.Header,'CreationDate')
truet = Expt.Header.CreationDate + t ./(24 .* 60 .* 60 .* 10000);
else
truet = Expt.Header.builddate + t ./(24 .* 60 .* 60 .* 10000);
end

function d = FindElectrodeDepth(Expt, t)

bid = find(Expt.Trials.Start < t);
aid = find(Expt.Trials.Start > t);
if isempty(aid) && isempty(bid)
    id = 1;
elseif isempty(aid)
    id = bid(end);
elseif isempty(bid)
    id = aid(1);
else
    id = bid(end);
end
d = Expt.Trials.ed(id);

function nt = CountTrials(Expt, t)
bid = find(Expt.Trials.Start > t(1) & Expt.Trials.Start < t(2));
nt = sum(Expt.Trials.Result(bid) > 0);


