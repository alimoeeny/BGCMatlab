function [same, score] = SameCell(xc, varargin)
%same = SameCell(xc, varargin)
%From efficacy calculation, determine if the two cells are the same
%0 = not
%1 = true
%2 suggests that one cell is a subset of the other.
%see also CalcEfficacy
same = 0;

if length(xc.efficacy) > 5
    effid = [5 6]; %adjusted based on rate in central ms of xcorr
else
    effid = [1 2];
end
a = max(xc.efficacy(effid));
c = min(xc.efficacy(effic));
b = max(xc.efficacy([3 4])); %next peak. Should be lower

if isfield(xc,'p') && xc.p(1) == xc.p(2)
    score = NaN;
    same = 1;
elseif a > 0.1 && length(xc.efficacy) > 7 %corrected for nearby rate etc
%if eff high and higher than rate for whoel central ms its good
%Unlesss there is just one co-incidnece in the ms (eff(8)). 
%of if max co-incidennce in other bins is similar. eff(7) is ratio
%sum(central 3 bins) /max(other bins). If eff(8) is low, this can be too.
%use criterion that uses eff(8) also instead of 3? 
     if (a > b * 1.5 || a > 0.8) && X.efficacy(7) > 3 && X.efficacy(8) > 2
         same = 1;
         if a > 2 * c
             same = 2;
         end
     else
         fprintf('Dubious Efficacy %.1f<->%.1f, %.3f - also get %.3f nearby. Ratio for coincidence %.1f\n',...
             xc.probe(1),xc.probe(2),a,b,xc.efficacy(7));
     end
elseif a > 0.1 && a > b * 1.5 %should only happen with old efficacy calculations
    same = 1;
else
end

