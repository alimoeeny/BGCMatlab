function [mev, details] = PlotExptEM(Expt, varargin)
%
% PlotExptEM plots eye movement data (and loads the data into Expt if
% needed)
%
% PlotExptEM(Expt, 'choice') returns mean eye positions by choice
% PlotExptEM(Expt, 'choicediff') returns difference in mean conjugate eye
% positions by choice (DOWN/RIGHT - UP/LEFT)
% PlotExptEM(Expt, 'choicediff',pidx,nidx) gives ids for trials to calculate difference for
% (difference is nidx - pidx)
% In both cases, there is no adjustment for the neurons prefereed.   
%
trid = 1:length(Expt.Trials);
showplot = 1;
startzero = 0;
details.n = 0;
splitchoice = 0;
plotdiff = 0;
duration = 20000;
emskip = 0;
j = 1;
while j <= nargin-1
    str = varargin{j};
    if strncmpi(str,'choice',3)
        splitchoice = 1;
        if length(varargin) > j+1 & isnumeric(varargin{j+1})
            j = j+1;
            pid = varargin{j};
            j = j+1;
            nid = varargin{j};
        else
            pid = find([Expt.Trials.RespDir] < 0);  %UP/LEFT
            nid = find([Expt.Trials.RespDir] > 0); %DOWN/RIGHT
        end
        if strfind(str,'diff')
            plotdiff = 1;
        end
    elseif strncmpi(str,'duration',3)
        j = j+1;
        duration = varargin{j};
    elseif strncmpi(str,'emskip',3)
        j = j+1;
        emskip = varargin{j};
    elseif strncmpi(str,'trials',3)
        j = j+1;
        trid = varargin{j};
    elseif strncmpi(str,'zero',3)
        startzero = 1;
    end
    j = j+1;
end

if ~isfield(Expt.Trials,'EyeData')
    Expt = LoadEMData(Expt);
end
sample_rate = 1/(10 * Expt.Header.CRsamplerate);
eyevals = cat(3,Expt.Trials(trid).EyeData);
mev = mean(eyevals,3);
sds = std(mev);
if sds(2) < 1e-10 & sds(4) < 1e-10
    details.eyes = [1 0 1 0];
else
    details.eyes = [1 1 1 1];
end

times = ([1:size(mev,1)] .* 0.1./Expt.Header.CRsamplerate); %in ms
if isfield(Expt.Header,'emtimes')
times = times + Expt.Header.emtimes(1)./10;
times = Expt.Header.emtimes./10;
end
details.times = times;

if splitchoice & isfield(Expt.Trials,'Saccades')
 pm = PlotExptEM(Expt,'Trials',pid);
 psacs = [Expt.Trials(pid).Saccades];
 nsacs = [Expt.Trials(nid).Saccades];
 id = find([psacs.start] < duration & [psacs.start] > emskip);
 [x,y] = pol2cart([psacs(id).dir], [psacs(id).size]);
 details.prefsac = mean(x) + i * mean(y);
 details.prefsacs = x + i * y;
 details.prefsact = [psacs(id).start];
 id = find([nsacs.start] < duration & [nsacs.start] > emskip);
 [x,y] = pol2cart([nsacs(id).dir], [nsacs(id).size]);
 details.nullsac = mean(x) + i * mean(y);
 details.nullsacs = x + i * y;
 details.nullsact = [nsacs(id).start];
 nm = PlotExptEM(Expt,'Trials',nid);
 %remove from each choice the overall starting mean ACROSS all choices.
   sid = find(times > 100 & times < 300);
    for j = 1:size(mev,2)
        pm(:,j) = pm(:,j) - mean(mev(sid,j));
        nm(:,j) = nm(:,j) - mean(mev(sid,j));
    end
    mev = cat(3,pm, nm);
    if plotdiff
        mev = diff(mev,1,3);
        if details.eyes(2) == 0 %r eye bad
            emv(:,1) = mev(:,1);
            emv(:,2) = mev(:,3);
        else
            emv(:,1) = mean(mev(:,[1 2]),2);
            emv(:,2) = mean(mev(:,[3 4]),2);
        end
        mev = emv;
    end
else 
    pm = [];
end


if startzero
    for j = 1:size(mev,2)
        mev(:,j) = mev(:,j) - mean(mev(1:10,j));
    end
end

if showplot
    if splitchoice & length(pm)
        if plotdiff
            plot(times,pm-nm);
        else
            plot(times,pm);
                hold on;
            plot(times,nm,':');
        end
    else
        plot(times,mev);
        legend('LH','RH','LV','RV');
    end
end