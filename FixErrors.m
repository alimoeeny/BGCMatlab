function [Errors, fixed] = FixErrors(name, varargin)
%FixErrors(dir|name) Read Error file and tidy up old errors

Errors = [];
fixed = 0;
savefixes = 1;

if ischar(name)
    if isdir(name)
        d = name;
        name = [d '/Errors.mat'];
        Errors = FixErrors(name);
        return;
    elseif exist(name)
        load(name);
        [Errors, fixed] = FixErrors(Errors, varargin{:});
        if fixed && savefixes
            fprintf('Saving Changes to %s\n',name);
            save(name,'Errors');
        end
        fullv.ShowErrors(Errors,'-all');
        return;
    end
elseif isstruct(name) && isfield(name,'s')
    Errors = name;
    fixed = 0;
    for j = 1:length(Errors);
        prog = '';
        type = '';
        eid = 0;
        severity = NaN;
        if isfield(Errors,'program') && ~isempty(Errors(j).program)
            prog = Errors(j).program;
        end
        if isfield(Errors,'progname') && ~isempty(Errors(j).progname)
            prog = Errors(j).progname;
        end
        exptno = 0;
        if isfield(Errors,'eid') && ~isempty(Errors(j).eid)
            exptno = Errors(j).eid;
        end
        if isfield(Errors,'expt') && ~isempty(Errors(j).expt)
            exptno = Errors(j).expt;
        end
        if exptno > 0 && (~isfield(Errors,'exptno') || isempty(Errors(j).exptno))
            Errors(j).exptno = exptno;
        end
        tprog = '';
        if strfind(Errors(j).s,'StimOns before')
            %%Usually just means ran a bunch of trials at start before
            %%storings
            tprog = 'APlaySpkFile';
            severity = 0;
            type = 'storage';
        elseif strncmp(Errors(j).s,'Missing StartExpt',15)
            % Usually just means that startexpt code was sent by binoc
            %too early for SPike2 to save it.  Should be gone after???
            tprog = 'APlaySpkFile';
            severity = 0;
            type = 'storage';
        elseif regexp(Errors(j).s,' [0-9]+ Chans at .* Excluded Time Range')
            tprog = 'BuildFullV';
            severity =1; %not necessarily in Trial
            type='datamissing';
        elseif regexp(Errors(j).s,'Missing [0-9]+ Trials for')
            tprog = 'LoadFullV';
            severity = 1;
            type = 'trialstructure';
        elseif regexp(Errors(j).s,'[0-9]+ Frames in.*bnc')
            tprog = 'APlaySpkFile';
            severity = 1;
            type = 'trialstructure';
        elseif regexp(Errors(j).s,'Trial [0-9]+.*missing probes')
            tprog = 'BuildFullV';
            severity =2;
            type='datamissing';
        elseif regexp(Errors(j).s,'Trial [0-9]+.*missing probes')
            tprog = 'BuildFullV';
            severity =2;
            type='datamissing';
        elseif regexp(Errors(j).s,'Duration varies')
            tprog = 'expt.Check';
            severity =1;
            type='duration';
        elseif regexp(Errors(j).s,'blocks end \(n\)')
            tprog = 'BuildFullV';
            severity =1;
            type='spkblk';
        elseif strncmp(Errors(j).s,'Two Frame counts',15)
            tprog = 'SetExptRC';
            severity =1;
            type='RCframecount';
        end
        if isempty(prog)
            prog = tprog;
        end
%This defines "Canonical" form for errs        
        if (~isfield(Errors,'progname') || isempty(Errors(j).progname)) && ~isempty(prog)
            Errors(j).progname = prog;
            fixed = fixed+1;
        end
        if (~isfield(Errors,'type') || isempty(Errors(j).type)) && ~isempty(type)
            Errors(j).type = type;
            fixed = fixed+1;
        end
        if (~isfield(Errors,'severity') | isnan(Errors(j).severity)) & ~isnan(severity)
            Errors(j).severity = severity;
            fixed = fixed+1;
        end
    end
    [strs, uid] = unique({Errors.s});
    if length(uid) < length(Errors)
        Errors = Errors(uid);
        fixed = fixed+1;
    end
    Errors = rmfields(Errors,{'program' 'eid' 'expt'});
elseif isstruct(name) && isfield(name,'errs') 
elseif isstruct(name) && isfield(name,'errmsg')
    name.Errors = name.errdata;
    for j = 1:length(name.errmsg)
        name.Errors(j).s = name.errmsg{j};
    end
    name.Errors = FixErrors(name.Errors);        
        
        
end