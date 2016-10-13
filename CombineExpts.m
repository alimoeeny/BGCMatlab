function Expt = CombineExpts(varargin)

j = 1;
nx = 0;
sorttrials = 0;
psychonly = 0;
while j <= length(varargin)
    if iscell(varargin{j}) && (isstruct(varargin{j}{1}) ||isstruct(varargin{j}{end}))
        Expts = varargin{j};
        nx = length(Expts);
    elseif iscell(varargin{j}) && sum(isexpt(varargin{1}{1})) %cell array of cell arrays
        dates = [];
        for k = 1:length(varargin{j})
            Expts{k} = CombineExpts(varargin{j}{k});
            if isfield(Expts{k}.Header,'CreationDate')
                dates(k) = Expts{k}.Header.CreationDate;
            end
        end
        if length(dates) == length(Expts)
            [~,id] = sort(dates);
            Expts = Expts(id);
        end
        nx = length(Expts);
    elseif isstruct(varargin{j}) & isfield(varargin{j},'Trials');
        nx = nx+1;
        Expts{nx} = varargin{j};
    elseif strcmp(varargin{j},'sort');
        sorttrials = 1;
    elseif strncmp(varargin{j},'psych',5) %psych only
        psychonly = 1;
    end
    j = j+1;
end

goodx = [];
nt = 1;
for j = 1:nx
    if isfield(Expts{j},'Trials') 
        goodx = [goodx j];
        if isfield(Expts{j}.Trials,'id') 
            if ~isempty(Expts{j}.Trials(1).id)
                fid(j) = Expts{j}.Trials(1).id;
            else
                fid(j) = min([Expts{j}.Trials.id]);
            end
        else
            fid(j)  = nt;
        end
        nt = nt + length(Expts{j}.Trials);
    end
end
nx = length(goodx);
Expts = Expts(goodx);
fid = fid(goodx);
if sorttrials
    [~,id] = sort(fid);
    Expts = Expts(id);
end

%For Stimvals that change between blocks, need to add these to each tril
%but only if its a meaninful trial variable
for k = 1:nx
    if ~isfield(Expts{k}.Trials,'Trial')
        fprintf('Filliing Trials.Trial in Expt%d\n',k);
        for j = 1:length(Expts{k}.Trials)
            Expts{k}.Trials(j).Trial = j;
        end
    end
end
fn = fields(Expts{1}.Stimvals);
fn = setdiff(fn,{'nskipped' 'i3' 'Pd' 'nt' 'rw'});
diffs = zeros(size(fn));
for j = 1:length(fn)
    for k = 2:nx
        if isnumeric(Expts{1}.Stimvals.(fn{j})) & isfield(Expts{k}.Stimvals,fn{j})
            if Expts{k}.Stimvals.(fn{j}) ~= Expts{1}.Stimvals.(fn{j})
                diffs(j) = 1;
            end
        end
    end
end
id = find(diffs);
for j = 1:length(id)
    filled(j) = 0;
    for k = 1:nx
        if ~isfield(Expts{k}.Trials,fn{id(j)})
            if ~filled(j)
                fprintf('Filling %s Ex%d\n',fn{id(j)},k);
                filled(j) = 1;
            end
            Expts{k} = FillTrials(Expts{k},fn{id(j)});
        end
    end
end


if ~isfield(Expts{1}.Header,'BlockStart')
    Expts{1}.Header.BlockStart = Expts{1}.Trials(1).Trial;
end
Expt = Expts{1};
for x = 2:nx;
    
    nt = Expt.Trials(end).Trial;
    ni = length(Expt.Trials);
    if Expts{x}.Trials(1).Trial < nt
        no = Expts{x}.Trials(1).Trial - nt +1;
    else
        no = 0;
    end
    for j = 1:length(Expts{x}.Trials)
        Expts{x}.Trials(j).Trial = Expts{x}.Trials(j).Trial + no;
        f = fields(Expts{x}.Trials);
        for k = 1:length(f)
        Expt.Trials(ni+j).(f{k}) = Expts{x}.Trials(j).(f{k});
        end
    end
%    Expt.Trials = [Expt.Trials Expts{x}.Trials];
     if isfield(Expts{x}.Header,'BlockStart')
         Expt.Header.BlockStart = [Expt.Header.BlockStart Expts{x}.Header.BlockStart+nt];
     else
         Expt.Header.BlockStart(end+1) = Expt.Trials(ni+1).Trial;
     end
end
Expt.Header.Combineids = GetExptNumber(Expts);
Expt.Header.suffixes = GetExptNumber(Expts);
f = fields(Expt.Trials);
if psychonly && isfield(Expt.Trials,'RespDir')
    for j = 1:length(Expt.Trials)
        if isempty(Expt.Trials(j).RespDir)
            Expt.Trials(j).RespDir = 0;
        end
    end
    id = find(ismember([Expt.Trials.RespDir],[-1 1]));
    Expt.Trials = Expt.Trials(id);
end
    
%Expt.Header.BlockStart(end) = length(Expt.Trials);


      