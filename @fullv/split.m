function split(FullV, varargin)
%fullv.split(FullV, 'trial', N) Creates two FullV files from one
%and writes out ExptXXaFullV.mat for the second half. The second
%file start with trial N. N.B. this is just the Nth trial in the Expt,
%not a trial id.

Expt = [];
splittrial = 0;
overwrite = 1;
excludeid = [];

if ischar(FullV)
    filename = FullV;
    FullV = LoadFullV(filename,'noconvert');
end

tsplit = 0;
eid = GetExptNumber(FullV);

j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'Header')
        Expt = varargin{j};
    elseif strncmpi(varargin{j},'trial',4)
        j = j+1;
        splittrial = varargin{j};
    end
    j = j+1;
end


if splittrial
    if isempty(Expt)
        Expts = ReadExptDir(fileparts(FullV.loadname));
        Expt = Expts{eid};
    end
    t(1) = Expt.Trials(splittrial-1).End(end);
    t(2) = Expt.Trials(splittrial).Start(1);
    tsplit = mean(t)./10000;
end

if tsplit > 0
    X = FullV;
    [a,b] = min(abs(FullV.t - tsplit));
    FullV.t = FullV.t(b:end);
    FullV.V = FullV.V(:,b:end);
    FullV.meanV = FullV.meanV(b:end);
    FullV.firstblk = 1;
    FullV.exptno = FullV.exptno+0.1;
end

outname = strrep(FullV.loadname,'FullV','aFullV');
if ~exist(outname) || overwrite
    fprintf('Saving %s\n',outname);
    save(outname,'FullV');
    suffs = {'ClusterTimes.mat' 'ClusterTimesDetails.mat' 'AutoClusterTimes.mat' 'AutoClusterTimesDetails.mat'}
    for j = 1:length(suffs)
        cfile = strrep(outname,'FullV.mat',suffs{j});
        if ~exist(cfile)
            afile = strrep(outname,'aFullV.mat',suffs{j});
            if exist(afile)
                fprintf('Copying %s to %s. Need to reclassify spikes\n',afile,cfile);
            end
        else
            fprintf('%s Already exists. Check that its been reclassified\n',cfile);
        end
    end
end
FullV = X;
FullV.lastsample = b;
fprintf('Saving %s\n',FullV.loadname);
save(FullV.loadname,'FullV');

