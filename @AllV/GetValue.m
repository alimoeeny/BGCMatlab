function value = GetValue(DATA, type, varargin)
%value = AllV.GetValue(DATA, type, varargin)
%  'layoutfile', 'last'  builds path to last saved layout
%  'configfile'  'last' builds path to last saved config

value = [];
strargs = cell2cellstr(varargin);

if strcmp(type,  'layoutfile')
    if isempty(DATA.layoutfile) || ~exist(DATA.layoutfile)
        value = DATA.defaultlayout;
        if DATA.options.UseLastLayout || sum(strcmp('last',strargs));
            value = strrep(DATA.defaultlayout,'.layout','last.layout');
        end
    else
        value = DATA.layoutfile;
    end
elseif strcmp(type, 'fitnumber')
    if isfield(DATA,'cluster') && isfield(DATA.cluster,'autofits')
        value = AllV.GetValue(DATA.cluster, type, varargin{:});
    elseif isfield(DATA,'autofiti')
        if DATA.autofiti(1) ==0 && isfield(DATA,'bestfit')
            if length(DATA.bestfit) ==1
                value = DATA.bestfit;
            else
                value = floor(DATA.bestfit(1)) + DATA.bestfit(2)./10;
            end
        elseif length(DATA.autofiti) > 1
            value = floor(DATA.autofiti(1)) + DATA.autofiti(2)./10;
        else
            value = DATA.autofiti(1);
        end
    end
elseif strcmp(type, 'configfile')
    if sum(strcmp('last',strargs))
        value = strrep(DATA.defaultconfig,'.config','last.config');
    elseif isempty(DATA.configfile)
        value = DATA.defaultconfig;
    else
        value = DATA.configfile;
    end  
end