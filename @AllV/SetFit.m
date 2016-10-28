function fit = SetFit(fit, type, varargin)
%fit AllV.SetFit(fit, type) modifies the fit struct to reflect alternate
%fits...
%AllV.SetFit(fit,'slim') uses clusters that have been combined if
% there were too clse
%if type is numeric, retutns the nth alternate fit (n = type)
%see also AllV.CombineFits

if iscell(fit) 
    fits = fit;
    if strcmp(type,'consensus');
        for j = 1:length(fits)
            if isfield(fits{j},'alternates')
                for k = 1:length(fits{j}.alternates)
                    if strcmp(fits{j}.alternates{k}.type,'consensus')
                        fit = AllV.SetFit(fits{j}, 'consensus');
                        return;
                    end
                end
            end            
        end
    end
    if isnumeric(type) && ~isempty(type)
        fit = fits{floor(type(1))};
        if length(type) ==2
            fit = AllV.SetFit(fit,type(2));
        elseif rem(type,1) > 0
            fit = AllV.SetFit(fit,round(rem(type(1),1).*10));
        end
        return;
    end
end
if ~isfield(fit,'alternates')
    return;
end

if isnumeric(type) && length(fit.alternates) >= type
    if type == 0;
        return;
    end
    type = fit.alternates{type}.type;
end

for j = 1:length(fit.alternates)
    if isfield(fit.alternates{j},'type') && strcmp(fit.alternates{j}.type,type)
        fit = FixFit(fit, fit.alternates{j});
        fit.alternateid = j;
        fit.alternatetype = type;
        fit.fitnumber = fit.fitnumber+0.1 .* fit.alternateid;
    end
end




function sfit = FixFit(fit, fix)
sfit = fit;
if isfield(fix,'substitutions')
    cls = unique(fit.clst);
    id = find(ismember(fix.substitutions(:,1),fix.substitutions(:,2)));
    [~,sid] = sort(sum(ismember(fix.substitutions(:,1),fix.substitutions(:,2)),2));
   for j = 1:size(fix.substitutions,1) 
       a = fix.substitutions(sid(j),:);
       if a(2) > 0
           sfit.clst(sfit.clst == a(1)) = a(2);
       elseif a(2) < 0
%a(2) = -1 (-2) means move cluster down one (2), not necessarily combining.
%Put at end of list in FixFit, so just do in order here, in case a later
%modification does another substitution. (e.g. a manual combinecls)
           sfit.clst(sfit.clst == a(1)) = a(1)+a(2);
           sfit.space(a(1)+a(2),:) = sfit.space(a(1),:);
           sfit.space(:,a(1)+a(2)) = sfit.space(:,a(1));
       end
   end
   for j = 1:size(fix.substitutions,1) 
       a = fix.substitutions(sid(j),:);
       if a(2) < 0
       end
       
   end
   cls = unique(fit.clst);
   sls = unique(sfit.clst);
   id = find(ismember(cls,sls));
   sfit.D = MatrixPermute(fit.D,id);
   if isfield(fit,'nearD')
       sfit.nearD = MatrixPermute(fit.nearD,id);
   end
   needf = {'isolation' 'NDisolation'};
   for j = 1:length(needf)
       f = needf{j};
       if isfield(fix,f)
           sfit.(f) = fix.(f);
       elseif isfield(fit,f)
           sfit.(f) = fit.(f)(id,:);
       end
   end
   needf = {'SU' 'dropi' 'bestscores' 'energy' 'lratio' 'optellipse' 'ellipse' 'isolation_distance'};
   for j = 1:length(needf)
       f = needf{j};
       if isfield(fit,f)
           sfit.(f) = fit.(f)(id);
       end
   end
elseif isfield(fix,'clst')
   sfit.clst = fix.clst;
   cls = unique(fit.clst);
   if isfield(fix,'isolation')
       sfit.isolation = fix.isolation;
   end
   if isfield(fix,'nocluster')
       sls = unique(sfit.clst);
       nid = find(sls == fix.nocluster);
       sfit.clst(sfit.clst == fix.nocluster) = 1;
   else
       sls = cls;
   end
%if new clst hase eliminated some clusters, adjust measures too
   [~, cid] = intersect(cls,sls);
   sfit.cls = cls(cid);
   needf = {'isolation' 'NDisolation'};
   for j = 1:length(needf)
       f = needf{j};
       if isfield(fix,f)
           sfit.(f) = fix.(f);
       elseif isfield(fit,f)
           sfit.(f) = fit.(f);
       end
       sfit.(f) = sfit.(f)(cid,:);
   end
%       sfit.isolation(nid,:) = [];
%       sfit.SU(nid) = [];
   needf = {'SU' 'dropi' 'bestscores' 'energy' 'lratio' 'optellipse' 'ellipse' 'isolation_distance'};
   for j = 1:length(needf)
       f = needf{j};
       if isfield(fit,f)
           sfit.(f) = fit.(f)(cid);
       end
   end
   sfit.D = MatrixPermute(fit.D,cid);
   if isfield(fit,'nearD')
       sfit.nearD = MatrixPermute(fit.nearD,cid);
   end

end

%better to keep it surely, with ellipses left matching original clusters? 
%   sfit = rmfields(sfit,'ellipse');
   if isfield(fix,'ellipse')
       sfit.ellipse = fix.ellipse;
   end
