function [DATA, details]  = JamesAutoCut(DATA, varargin)usev = 0;retrigger = 0;reapply = 0;verbose = 0;j = 1;while j <= length(varargin)    if strncmpi(varargin{j},'reapply',5)        reapply= 1;        j = j+1;        C = varargin{j};    elseif strncmpi(varargin{j},'retrigger',5)        retrigger = 1;    end    j = j+1;endaddpath('/b/bgc/matlab/james/autocluster');Vall = AllV.mygetappdata(DATA,'Vall');Vall.Vtime = Vall.t;Vall.Fs = 40000;X.Vtime = Vall.t;X.Fs = 40000;X.V = double(Vall.V').*Vall.intscale(1)/Vall.intscale(2);clust_params.gmm_inits = 100;clust_params.min_Pcomp = 0.05;clust_params.use_ms = 0;clust_params.try_features = [1 2 4];clust_params.max_back_comps = 2;clust_params.cluster_biase = 0.5;clust_params.target_rate = 50;if DATA.interactive < 0    clust_params.summary_plot = 0;else    clust_params.summary_plot = 1;end  if retrigger    [details,spike_features,sum_fig] = detect_and_cluster_init(X,clust_params,DATA.chspk);    uids = setdiff(1:length(details.spk_inds),details.artifact_ids);    details.spike_xy = details.spike_xy(uids,:);    details.spike_clusts = details.spike_clusts(uids);    details.times = details.times(uids); %spk inds is done below because it also needs doing after reapply %(?because?)    DATA.uid = 1:length(uids);    spike_features = spike_features(uids,:);    elseif reapply    [details,spike_features,spike_xy,Spikes] = apply_clustering(X, C.jamescluster, C.jamescluster.params, 0, []);%needs to return     [cl, cluster, xy]    DATA.uid = 1:length(details.spike_clusts);    DATA.clst = details.spike_clusts;    DATA.xy{1} = spike_xy;    uids = 1:length(details.spk_inds); %artifacts already outelse    Spikes = AllV.MakeJamesSpikes(DATA);    [details,spike_features,sum_fig] = detect_and_cluster_init(X,clust_params,DATA.chspk,Spikes);    DATA.uid = details.uids;    DATA.clst = zeros(size(DATA.clst));    spike_features(details.uids,:) = spike_features;    spike_features(details.artifact_ids,:) = NaN;    details.spike_xy(details.uids,:) = details.spike_xy;    details.spike_xy(details.artifact_ids,:) = NaN;endDATA.cluster.autocutmode = 'james';DATA.cluster.cluster = 1;C = DATA.cluster;gid = 1:length(details.spk_inds);if reapply || retrigger    Tv = AllV.mygetappdata(DATA,'TriggerV');    Vall = AllV.mygetappdata(DATA,'Vall');    vid = details.spk_inds(uids);    DATA.xy{1} = details.spike_xy;    if DATA.usetrials        [gid, postid] = TimesInTrial(Vall.t(vid), DATA.Expt, DATA.preperiod, DATA.postperiod);        vid = vid(gid);        DATA.xy{1} = details.spike_xy(gid,:);    end    DATA.uid = 1:length(vid);    DATA.rV = Tv(vid);    DATA.t = Vall.t(vid);    details.comp_idx = details.comp_idx(gid);    DATA.jamescluster = rmfields(details,{'spk_inds' 'uids' 'times' 'spike_clusts' 'spike_xy'});    DATA.cluster.jamescluster = DATA.jamescluster;    DATA.cluster.ctime = now;    [AllVoltages, DATA] = AllV.BuildAllV(DATA, vid', details.params.spk_pts);    DATA.spikefeatures = spike_features(gid,:);    DATA = AllV.SetPCs(DATA,0,0);    DATA.clst = details.spike_clusts(gid);    cls = unique(DATA.clst);    fprintf('Final James Cut %d clusters (including hash)\n',sum(cls>0));    DATA.spikefeatures = cat(2,spike_features(gid,:),details.spike_xy(gid,:));    DATA.clid = find(DATA.clst == C.cluster+1);    DATA.nid = find(DATA.clst ~= C.cluster+1);    X = AllV.CalcMyDistanceMatrix(DATA,DATA.cluster,verbose,'nearmahal');    DATA.cluster.autofits{1} = X;    DATA.cluster.jamescluster = CopyFields(DATA.cluster.jamescluster, X, {'isolation' 'D' 'nearD'});endDATA.gmcid = details.comp_idx;C = DATA.cluster;C.next = {};C.probe = DATA.probe;C.space = [13 1 2];C.shape = 3;C.chspk = DATA.chspk;C.auto  = 2;C.angle = 0;C.autocutmode = DATA.autocutmode;C.Trigger = details.trig_thresh;C.MeanSpike = AllV.PlotMeanSpike(DATA,'recalc');C.bmc = 0;C.mahal(1) = details.dprime;C.mahal(2) = details.iso_dists(1);C.mahal(3) = details.Lratios(1);C.mahal(4) = details.Lratios(1);C.ctime = now;for j = 2:length(details.Lratios)    C.next{j-1}.mahal(3) = details.Lratios(j);    C.next{j-1}.mahal(2) = details.iso_dists(j);    C.next{j-1}.mahal(1) = details.dprime;enddetails.newDATA = 0;details.clust_params = clust_params;C.cluster = 1;C.ncut = length(DATA.clid);C.nspks = length(DATA.clst);C = rmfields(C,'eckercluster'); %remove any old C.autocutmode = 'james';DATA.cluster = C;