function Expt = FixEd(Expt)


id = find(Expt.Trials.ed ~= 0)
if id(1) > 10 & diff(Expt.Trials.ed(id(1)-1:id(1))) > 0.1
    Expt.Trials.ed(1:id(1)) = Expt.Trials.ed(id(1));
end
