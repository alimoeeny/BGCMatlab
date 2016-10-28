function I = BuildLFPIndex(name,varargin)
%I = BuildLFPIndex(name,varargin)
%For Utah Folders, read in all the lfp files that are made and summarize contents
%For Spike2 Folders, compare Expts List with list of .lfp files
listonly = 0;
j = 1;
defaultarray = '';

while j <= length(varargin)
    if strncmpi(varargin{j},'fixlist',4)
        listonly = 2;
    elseif strncmpi(varargin{j},'list',4)
        listonly = 1;
    end
    j = j+1;
end

if iscellstr(name)
    for j = 1:length(name)
        I{j} = BuildLFPIndex(name{j},varargin(:));
    end
    return;
elseif iscell(name) %previous result
    for j = 1:length(name)
        if isfield(name{j},'D')
            D = name{j}.D;
            nerr = 0;
            errs = {};
            for k = 1:length(D)
                x = PrintErrorSummary(D{k});
                errs = {errs{:} x{:}};
                nerr = nerr + length(x);
            end
            if nerr == 0
                fprintf('No errors in %s\n',name{j}.name);
            else
                PrintErrorSummary(errs,'show');
            end
        end
    end
    return;    
end
if strcmp(name,'jbe')  %buidl index for jbe Utah data
    d = mydir('/b/data/jbe/G*');
    I = BuildLFPIndex({d.name},varargin{:});
    return;
end

if isdir(name)
    mnk = GetMonkeyName(name);
    if isempty(defaultarray)
        if strcmp(mnk,'jbe')
            defaultarray = 'Utah';
        else
            defaultarray = 'uprobe';
        end
    end
    goodexpts = 0;
    outname = [name '/LFPlist.mat'];
    Array = GetArrayConfig(name);
    if isempty(Array) || ~isfield(Array,'type')
            Array.type = defaultarray;
    end
    if listonly 
        if sum(strncmpi(Array.type,{'uprobe' 'vprobe'},4))
            Expts = ReadExptDir(name);
            for j = 1:length(Expts)
                lfpname = [name '/' GetName(Expts{j},'lfp')];
                if ~exist(name)
                    fprintf('Missing %s\n',name);
                    goodexpts(j) = 0;
                else
                    goodexpts(j) = 1;
                end
            end
            I.name = outname;
            I.havefiles = goodexpts;
            I.missing = sum(goodexpts ==0);
            I.ngood = sum(goodexpts > 0);
        elseif strncmpi(Array.type,'utah',4)
            I.name = outname;
            if ~exist(outname)
                I.list = 0;
                fprintf('No LFP List for %s\n',name);
            else
                I.list = 1;
                I.exptno = [];
                I.lowpasserror = 0;
                OUT = load(outname);
                if ~iscell(OUT.D)
                    D{1} = OUT.D;
                    OUT.D = D;
                end
                for j = 1:length(OUT.D)
                    exptno = 0;
                    ncut = 0;
                    
                    if isfield(OUT.D{j},'exptno')
                        name = [outname 'Expt' num2str(OUT.D{j}.exptno)];
                        exptno = OUT.D{j}.exptno;
                    end
                    I.exptno = [I.exptno exptno];
                    if isfield(OUT.D{j},'lfpnames')
                        I.names{j} = OUT.D{j}.lfpnames{1};
                        name = OUT.D{j}.lfpnames{1};
                    elseif exptno > 0
                        name = [outname 'Expt' num2str(exptno)];
                    else
                        name = outname;
                    end
                    if isfield(OUT.D{j},'errs')
                        I.nerrs(j) = length(OUT.D{j}.errs);
                        for e = 1:length(OUT.D{j}.errs)
                            err = OUT.D{j}.errs{e};
                            if strfind(err,'points, array is')
                                fprintf('Bad Channel in %s\n',name);
                            elseif strfind(err,'filtered out')
                                ncut = ncut+1;
                            else
                                fprintf('Unknown Errors in Expt %d\n',j);
                            end
                        end
                    else
                        I.nerrs(j) = 0;
                    end
                    I.lowpasserror(j) = ncut;
                end
                if  sum(I.lowpasserror)
                    fprintf('Filteriing wrong in %s\n',name);
                end
            end
        else
            I.err = sprintf('Unknown Array Type %s',Arrat,type);
            cprintf('red','%s in %s\n',I.err,name);
            I.name = outname;
        end
        if listonly == 1
            return;
        end
    else
        I.name = outname;
    end
    if exist(outname)
        OUT = load(outname);
        if isfield(OUT,'D')
            if  isstruct(OUT.D)
                olddata = OUT.D;
            OUT.D = {};
            else
                D = {};
                for j = 1:length(OUT.D)
                    expts(j) = ExptIndex(OUT.D{j});
                    if expts(j) > 0
                        D(expts(j)) = OUT.D(j);
                    end
                end
                OUT.D = D;
            end
        end
        
    end
    if sum(strncmpi(Array.type,{'uprobe' 'vprobe'},4))
        Expts = ReadExptDir(name);
        I.LFPerrs = {};
        I.LFPerrdata =[];
        for j = 1:length(Expts)
            [~, E] = LoadLFP(Expts{j});
            if E.gotlfp == 0
                goodexpts(j) = 0;
            elseif isfield(E,'LFPerrs')
                goodexpts(j) = 2; %have data, but with errors
                I.LFPerrs = cat(1,I.LFPerrs(:),E.LFPerrs(:));
                I.LFPerrdata = cat(2,I.LFPerrdata,E.LFPerrdata);
            else
                goodexpts(j) = 1;
            end
        end
        I.goodexpts = goodexpts;
        if sum(goodexpts)
            save(outname,'-struct','I');
        end
    elseif strncmpi(Array.type,'Utah',4) 
        d = mydir([name '/Expt*.lfp.mat']);
        if ~isempty(d)
        for j = 1:length(d)
            X = load(d(j).name);
            e = GetExptNumber(d(j).name);
            if isfield(X,'LFPHeader') 
                OUT.D{e}.Header = X.LFPHeader;
                goodexpts(e) = 1;
                OUT.D{e}.npts = size(X.LFP.rawlfp,1);
                OUT.D{e}.nprobes = size(X.LFP.rawlfp,2);
            elseif isfield(X,'LFP')
                OUT.D{e}.npts = size(X.LFP.rawlfp,1);
                OUT.D{e}.nprobes = size(X.LFP.rawlfp,2);
                if isfield(X.LFP,'Header')
                    OUT.D{e}.Header = X.LFP.Header;
                else
                    OUT.D{e} = CopyFields(OUT.D{e},X.LFP,'kernel','sd','decimate','samper');
                end
                goodexpts(e) = 1;
            end
            OUT.D{e}.exptno = e;
            if OUT.D{e}.nprobes < length(Array.X)
                fprintf('Only %d probes in %s\n',OUT.D{e}.nprobes,d(j).name);
                goodexpts(e) = 2;
            end
        end
        I.D = OUT.D;
        if sum(goodexpts)
            bid = find(~goodexpts); %no data here any more
            for j =1:length(bid)
                OUT.D{bid(j)} = [];
            end
            save(outname,'-struct','OUT');
        end
        end
    else
        cprintf('red','Unknown Array Type %s in %s\n',Array.type,name);
    end
else
    I.name = name;
end


function LFPerrs = PrintErrorSummary(LFP,varargin)

show = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'show',4)
        show = 1;
    end
    j = j+1;
end

LFPerrs = {};
if isfield(LFP,'errs')
    for k = 1:length(LFP.errs)
        probes(k) = GetProbeFromName(LFP.errs{k});
        errs{k} = deblank(regexprep(LFP.errs{k},'.p[0-9]+FullV.mat',''));
    end
    [a,b] = unique(errs);
    for k = 1:length(a)
        id = find(strcmp(a{k},errs));
        LFPerrs{end+1}.s = a{k};
        LFPerrs{end}.probes = probes(id);
        if length(id) > 10
            pstr = sprintf('%d-%d',min(probes(id)),max(probes(id)));
        else
            pstr = sprintf(' %d',probes(id));
        end
        if show
            fprintf('%s Probes %s\n',a{k},pstr);
        end
    end
elseif isfield(LFP,'Header') && isfield(LFP.Header,'errs')
    LFPerrs = PrintErrorSummary(LFP.Header);
elseif iscell(LFP)
    for k = 1:length(LFP)
        expts(k) = GetExptNumber(LFP{k}.s);
        errs{k} = deblank(regexprep(LFP{k}.s,'Expt[0-9]+',''));
    end
    [a,b] = unique(errs);
    for k = 1:length(a)
            id = find(strcmp(a{k},errs));
            LFPerrs{end+1}.s = a{k};
            LFPerrs{end}.expts = expts(id);
        if show
            fprintf('%s Expts %s\n',a{k},sprintf(' %d',expts(id)));
        end
    end
        
    
end


function eid = ExptIndex(X)

eid = 0;
if isfield(X,'exptno')
    eid = X.exptno;
elseif isfield(X,'Header')
end
