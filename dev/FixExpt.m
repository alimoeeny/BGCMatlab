function Expt = FixExpt(Expt,type)


if strcmp(type,'ed') & isfield(Expt.Trials,'ed')

id = find(Expt.Trials.ed ~= 0);
if length(id) > 1 & id(1) > 10 & diff(Expt.Trials.ed(id(1)-1:id(1))) > 0.1
    Expt.Trials.ed(1:id(1)) = Expt.Trials.ed(id(1));
end
end
