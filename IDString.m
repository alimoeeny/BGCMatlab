function s = IDString(X, varargin)
strargs = cell2cellstr(varargin);

s = [];
type = [];
j = 1; 
while j <= length(varargin)
    if strncmp(varargin{j},'probe',4)
        type = probe;
    end
    j = j+1;
end

if isfield(X,'Header') && (isfield(X.Header,'cellnumber') || isfield(X.Header,'probe') || isfield(X.Header,'CreationDate')) 
    name = '';
    if isfield(X(1),'name')
        s = BuildFileName(X(1),'shortname');
        if ~isempty(s)
            name = [s ' '];
        end
    end
    P = X(1).Header; %in case its an array
    if ~isfield(P,'cellnumber')
        P.cellnumber = 0;
    end
    if ~isfield(P,'clusterid')
        P.clusterid = 0;
    end
    if isempty(type) || strcmp(type,'probe')
        if isfield(P,'probe')
            pstr = sprintf(' P%.0f',P.probe);
            if mod(P.probe,1) > 0.01
                pstr = sprintf('P%.1f',P.probe);
            end
        elseif isfield(X,'probe')
            pstr = sprintf(' P%.0f',X.probe);
        else
            pstr = '';
        end
        if P.clusterid > 0
            pstr = sprintf('%s.%d',pstr,P.clusterid);
        end
        if P.cellnumber > 0
            s = sprintf('%sCell %d%s',name,P.cellnumber, pstr);
        else
            s = sprintf('%sMU%s',name, pstr);
        end
    end
    if sum(strcmp('trialcount',strargs))
        nt = [];
        if isfield(X,'Trials')
            nt = length(X.Trials);
        elseif isfield(X,'Data') && isfield(X(1).Data,'Trials')
            nt = length(X(1).Data.Trials);
        end
        if ~isempty(nt)
            s = sprintf('%s %d Trials',s,nt);
        end
    end
end