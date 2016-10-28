 function [voffset, ylim] = SetVOffset(DATA, AllSpikes, e, setprobe) % voffset = SetVOffset(DATA, AllSpikes, e,  set offsets so that channels % don't overlap (too much) % ..., setprobe) 4th argument means set for this probe using xch  ylim = []; if nargin ==3     setprobe = 0; end if isfield(AllSpikes,'xchans')     setprobe = AllSpikes.probe; end if length(setprobe) > 1     probelist = setprobe; else %really want chspk here.  But Maybe make calls set it more ofter     probelist = 1:DATA.nprobes; end if isempty(AllSpikes)     voffset = [1:DATA.nprobes] .* 2; %? could use meanspiek - get values at startup elseif setprobe(1) > 0 &&  (isfield(AllSpikes,'xchans') || iscell(AllSpikes) && isfield(AllSpikes{e,setprobe(1)},'xchans'))     if isfield(AllSpikes,'xchans')         X = AllSpikes;     else         X = AllSpikes{e,setprobe(1)};     end     voffset(setprobe(1)) = diff(X.VRange);     if isfield(X,'xVranges') %have range for each extra probe         for j = 1:length(X.xchans)             if ismember(X.xchans(j),probelist)                 voffset(X.xchans(j)) = diff(X.xVranges(j,:));             end         end     else         for j = 1:length(X.xchans)             if ismember(X.xchans(j),probelist)                 voffset(X.xchans(j)) = diff(X.xVrange);             end         end     end     meanstep = diff(X.xVrange);     xid = setdiff(setprobe,[X.xchans X.probe]);     for j = 1:length(xid)         voffset(xid(j)) = meanstep;     end     voffset = cumsum(voffset);     vRange(:,1) = X.xVrange(1);     vRange(:,2) = X.xVrange(2);     if length(setprobe) > 1         allprobes = setprobe;         voffset = voffset - voffset(min(setprobe));     else         allprobes = [X.xchans setprobe(1)];         voffset = voffset - voffset(setprobe(1));     end     voffset = voffset - voffset(setprobe(1));          vmins = voffset+vRange(:,1)';     vmaxs = voffset+vRange(:,2)';     voffset(end+1:DATA.nprobes) = voffset(end);     voffset = voffset * 0.7;     ylim(1) = min(vmins(allprobes)) .* 0.7;     ylim(2) = max(vmaxs(allprobes)) .* 0.7;     ylim(1) = ylim(1) - diff(X.xVrange)*0.2; %avoid clipping at bottom elseif iscell(AllSpikes)     maxv = CellToMat(AllSpikes(e,:),'VRange');     xch = CellToMat(AllSpikes(e,:),'xchans');     if isempty(maxv)         voffset = [1:DATA.nprobes] .* 2; %? could use meanspiek - get values at startup     else         voffset =  cumsum(cat(1,maxv(2:end,2), maxv(end,2))) - cumsum(maxv(:,1));         if length(voffset) < max(xch(:))             X = AllSpikes{e,length(voffset)};             if isfield(X,'xVrange')                 voffset(max(xch(:))) = max(voffset) + diff(X.xVrange);             else                 voffset(max(xch(:))) = max(voffset) + 2.*X.xmaxv; %? 2* is too much?             end         end         voffset = voffset * 0.7;     end else     v = diff(AllSpikes.VRange);     if isempty(v)         v = 0;     end     voffset(setprobe) = v;     voffset(setprobe+1:DATA.nprobes) = voffset(setprobe); end