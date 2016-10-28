function DATA = SetData(DATA, fcn, varargin)
%DATA = AllV.SetData(DATA, ...)  make sure data has loaded relevant stuff
%       AllV.SetData(DATA, 'autofits') gets autofits from ClusterDetails

if strcmp(fcn, 'autofits')
   if ~isfield(DATA.cluster,'autofits')
        [DATA, CD] = AllV.LoadClusterDetails(DATA,'auto');
        if ~isempty(CD)
            DATA.cluster.autofits = CD{DATA.probe(1)}.autofits;
        end
    end  
elseif strcmp(fcn, 'checkellipsefits')
    for j = 1:length(DATA.cluster.autofits)
        if ~isfield(DATA.cluster.autofits{j},'ellipse')
            fit = DATA.cluster.autofits{j};
            [~,~, fit] = AllV.SetClusterFromFit(DATA, DATA.cluster, fit);
            if isfield(fit,'alternates')
                for a = 1:length(fit.alternates)
                    [~,~, fit] = AllV.SetClusterFromFit(DATA, DATA.cluster, fit, fit.alternates{a}.type);
                end
            end
            DATA.cluster.autofits{j} = fit;
        end
    end
end