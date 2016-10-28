function ShowErrors(res, varargin)


showall = 0;
prestr = [];
for j = 1:length(varargin)
    if strncmp(varargin{j},'-all',4)
        showall = 1;
    else
    prestr = [prestr varargin{j}];
    end
    j = j+1;
end

for j = 1:length(res)
    if isfield(res,'eid') && ~isempty(res(j).eid)
        estr = sprintf('E%d',res(j).eid);
    elseif isfield(res,'exptno') && ~isempty(res(j).exptno)
        estr = sprintf('E%d',res(j).exptno);
    elseif isfield(res,'expt') && ~isempty(res(j).expt)
        estr = sprintf('E%d',res(j).expt);
    else
        estr = '';
    end
    if isfield(res,'progname') && ~isempty(res(j).progname)
        progstr = [res(j).progname ': '];
    else
        progstr = '';
    end
    if isempty(progstr) || sum(strncmp(progstr,{'fullv'},4)) || showall
        if isfield(res,'errs')
            for e = 1:length(res(j).errs)
                fprintf('%s%s:%s%s\n',prestr,progstr,res(j).name,res(j).errs{e});
            end
        elseif isfield(res,'s')
            fprintf('%s%s%s: %s\n',progstr,prestr,estr,res(j).s);
        elseif isfield(res,'Errors')
            if isfield(res,'name')
                name = res(j).name;
            else
                name = '';
            end
            fullv.ShowErrors(res.Errors,name,varargin{:});
        elseif iscell(res{j})
            fullv.ShowErrors(res{j},varargin{:});
        elseif isfield(res{j},'errs')
            for e = 1:length(res{j}.errs)
                fprintf('%s:%s\n',res{j}.name,res{j}.errs{e});
            end
        end
    end
end