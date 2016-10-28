function [AllVoltages, details] = Spikes2V(S, Cluster, allnprobes, varargin)
%AllV.Spikes2V(S, Cluster)  Converts a Spikes file into AllVoltages Matrix
%
p = Cluster.probe(1);
if ~isfield(S,'probe')
    S.probe = Cluster.probe(1);
end
 if isfield(S.Header,'version')
        version = S.Header.version;
    else
        version =1.1;
    end

    if isfield(S,'xchans')
%Old files with n events less that n points transposed values, so xchans dimensiont
%didn't match
        if version < 1.14 && size(S.values,1) ~= length(S.times) %dont transpose
            AllVoltages = zeros([allnprobes size(S.values,1) size(S.values,2)]);
            AllVoltages(Cluster.probe(1),:,:) = double(S.values) .* S.maxv./S.maxint;
        else
            AllVoltages = zeros([allnprobes size(S.values,2) size(S.values,1)]);
            if isfield(S,'probe')
                AllVoltages(S.probe(1),:,:) = double(S.values') .* S.maxv./S.maxint;
            else
                AllVoltages(Cluster.probe(1),:,:) = double(S.values') .* S.maxv./S.maxint;
            end
        end
        chspk = Cluster.chspk;
        AllVoltages(S.xchans,:,:)= double(S.xvalues) .* S.xmaxv./S.maxint;
        newchspk = sort([S.probe S.xchans]);
        details.probes = newchspk;
        details.nprobes = length(newchspk);
        if sum(~ismember(chspk,newchspk)) %don't have channels
            details.missing = setdiff(chspk,newchspk);
            chspk = newchspk;
            fprintf('AllV.Spikes2V: Missing Spike Channels (%s) for%s Using %s\n',sprintf(' %d',details.missing),str(Cluster,'nocl'),str(S));
        end
    else
        AllVoltages(Cluster.probe(1),:,:) = double(S.values') .* S.maxv./S.maxint;
        chspk = Cluster.probe(1); %1 for single Fullv files
        details.nprobes = 1;
        details.probes = chspk;
    end
    details.t = S.times(:)'./10000; %column vector    
    details.chspk = chspk;