function [badtrial,badid, goodid] = FindMissingTrials(Expt, t, varargin)
%[badtrial,badid] = FindMissingTrials(Expt, t, varargin)
% find trials in Expt for which no elements in t fall in the trial

badtrial = [];
badid = [];
goodid = [];
trange = Expt.Trials(end).Start(1) - Expt.Trials(1).Start(1);
if diff(minmax(t(t>0))) < trange/1000 %wrong units
    t = t .* 10000;
end
tmax = max(t);
for j = 1:length(Expt.Trials)
    good(j) = sum(t> Expt.Trials(j).Start(1) & t < Expt.Trials(j).End(end));
    pregood(j) = sum(t> Expt.Trials(j).Start(1)-1000 & t <= Expt.Trials(j).Start(1));
    if good(j) == 0 || pregood(j) == 0
        badtrial(end+1) = j;
        badid(end+1) = Expt.Trials(j).id;
    elseif tmax < Expt.Trials(j).End(end)
        badtrial(end+1) = j;
        badid(end+1) = Expt.Trials(j).id;        
    else
        goodid(end+1,1) = find(t > Expt.Trials(j).Start(1),1);
        goodid(end,2) = find(t > Expt.Trials(j).End(end)-0.001,1);
    end
end
    