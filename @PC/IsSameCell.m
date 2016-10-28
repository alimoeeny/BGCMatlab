function [same, score] = IsSameCell(xc, varargin)
%[same, eff] = IsSameCell(xc, varargin)
%From efficacy calculation, determine if the two cells are the same
%0 = not
%1 = true
%2 suggests that one cell is a subset of the other.
%-1 suspicous # of coincidences, but does not look like real identity
%      can happen with very low counts
% eff returns the coincidences both ways, largest first
%see also CalcEfficacy
same = 0;
if isnumeric(xc) && length(xc) > 5
    a = xc;
    clear xc;
    xc.efficacy = a;
end
if length(xc.efficacy) > 5
    effid = [5 6]; %adjusted based on rate in central ms of xcorr
else
    effid = [1 2];
end
a = max(xc.efficacy(effid));
c = min(xc.efficacy(effid));
b = max(xc.efficacy([3 4])); %next peak. Should be lower
score = [a c];

if isfield(xc,'p') && xc.p(1) == xc.p(2)
    score = NaN;
    same = 1;
elseif a > 0.1 && length(xc.efficacy) > 7 %corrected for nearby rate etc
%if eff high and higher than rate for whoel central ms its good
%Unlesss there is just one co-incidnece in the ms (eff(8)). 
%of if max co-incidennce in other bins is similar. eff(7) is ratio
%sum(central 3 bins) /max(other bins). If eff(8) is low, this can be too.
%use criterion that uses eff(8) also instead of 3? 
     if (a > b * 1.5 || a > 0.8) && xc.efficacy(7) > 3 && xc.efficacy(8) > 2
         same = 1;
         if a > 2 * c
             same = 2;
         end
     else %fishy case
         same = -1;
     end
elseif a > 0.1 && a > b * 1.5 %should only happen with old efficacy calculations
    same = 1;
else
end

