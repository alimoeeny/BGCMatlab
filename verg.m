function verg(varargin)
%binoc
%GUI for running binoclean via pipe.

autoquit = 0;
uselastlayout = 0;
layoutfile = '';

if length(varargin) & ishandle(varargin{1})
    f = varargin{1};
    while ~isfigure(f)
        f = get(f,'parent');
    end
    DATA = get(f,'UserData');
    varargin = varargin(3:end);
else
    checkforrestart = 1;
    TOPTAG = 'vergwindow';
    it = findobj(allchild(0),'flat','Tag',TOPTAG,'type','figure');
    setdata = 0;
%just parse arguments here that need to be set before reading DATA from figure    
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'autoreopen',6)
        if j ==1 && ~isempty(it) %User changing after window opened
            DATA = get(it(1),'UserData');
            setdata = 1;
        end
        DATA.autoreopen =1;
        if length(varargin) == j+1          %final argument after autoreopen = line to resend
            j = j+1;
            DATA.reopenstr = varargin{j};
        end
    elseif strcmp(varargin{j},'autoquit')
        autoquit = 1;
    elseif strcmp(varargin{j},'host')
        j = j+1;
        DATA.ip = ['http://' varargin{j} ':1110/'];
    elseif strcmp(varargin{j},'lastlayout')
        uselastlayout = 1;
    elseif strcmp(varargin{j},'layout')
        j = j+1;
        layoutfile = varargin{j};
        if strcmp(layoutfile,'last')
            uselastlayout = 1;
        else
            d = dir('/local/*.layout');
            id = find(CellToMat(strfind({d.name},varargin{j})));
            if ~isempty(id)
                layoutfile = ['/local/' d(id(1)).name];
            end
        end
    elseif strcmp(varargin{j},'new')
        checkforrestart = 0;
    elseif strcmp(varargin{j},'record')
        DATA = OpenBinocLog(DATA,'both');
    elseif strcmp(varargin{j},'verbose')
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            DATA.verbose = varargin{j};
        elseif exist('DATA','var')
            DATA.verbose = ones(size(DATA.verbose));
        else
            verbose = 1;
        end
        DATA = OpenBinocLog(DATA,'frombinoc');
    elseif strcmp(varargin{j},'monitor')
        DATA = OpenBinocLog(DATA, 'frombinoc');
        DATA.perfmonitor = 1;
    end
    j = j+1;
end
if setdata
            set(it,'UserData',DATA);
end


it = findobj(allchild(0),'flat','Tag',TOPTAG,'type','figure');
if isempty(it)
    tt = TimeMark([], 'Start');
    DATA.tag.top = 'vergwindow';
    ts = now; 
    DATA.starttime = now;
    DATA = SetDefaults(DATA);
    vpath = mfilename('fullpath');
    DATA.name = strrep(DATA.vergversion,'verg.','Verg Ver ');
%open pipes to binoc 
% if a file is named in teh command line, then take all setting from there.
% Otherwise read from binoc
    vtest = {};
    if strcmp(varargin{1},'help')
        CodesPopup(DATA,[],'popup');
        return;
    elseif length(varargin) && exist(varargin{1},'file')
        DATA = OpenPipes(DATA, 0);
        if DATA.togglecodesreceived == 0
            cprintf('red','Binoclean did not send Toggle Codes within %.2f\n',DATA.readdur);
        end
        tt = TimeMark(tt, 'Pipes');
        DATA.stimfilename = varargin{1};
        PauseRead(DATA,1);
%why do we neeed to send state before reading stim files?
%.e. which parameters is this for. Could just send a short list...
        SendState(DATA); %params loaded from verg.setup, binoc.setup etc
        tt = TimeMark(tt, 'SendState');
        DATA = ReadStimFile(DATA, '/local/verg.setup','fromverg'); %make sure these go to binoc
        vtext = DATA.exptlines;
        tt = TimeMark(tt, 'VergSetup');
        if strcmp(varargin{1},'demo') 
            if isfield(DATA.verg,'demofile')
                varargin{1} = DATA.verg.demofile;
            else
                varargin{1} = '/local/c/binoclean/stims/demo/demo.stm';
            end
        end
                

        [DATA, details] = ReadStimFile(DATA,varargin{1}, 'init');
        DATA.state.stimfile = varargin{1};
        if details.badcodes > 1
            DATA.state.stimfileerrs = details.badcodes;
        end
        tt = TimeMark(tt, 'Read');
        SendCode(DATA,'vve'); %send verg version
        DATA.rptexpts = 0;  %can't set this in startup files, only quickmenus
        DATA = DrainBinocPipe(DATA);
        PauseRead(DATA,0);
        if exist(DATA.binoc{1}.lo,'file')
            DATA = ReadLogFile(DATA, DATA.binoc{1}.lo);
        end
        DATA = GetState(DATA,'stimfile');
        DATA = DrainBinocPipe(DATA);
        tt = TimeMark(tt, 'FromBinoc');

        j = 2;
        while j <= length(varargin)
            if ischar(varargin{j}) && sum(strncmp(varargin{j},{'new'},3)) == 0
                outprintf(DATA,'%s\n',varargin{j});
            end
            j = j+1;
        end
    else
        if ~isempty(varargin) && ischar(varargin{1})
            vergwarning(sprintf('Cant Read %s',varargin{1}));
        end
        DATA = ReadVergFile(DATA, DATA.layoutfile);
        DATA = OpenPipes(DATA, 1);
        DATA.nexpts = 1;
        DATA.Expts{1} = ExptSummary(DATA);
    end
    
    if uselastlayout
        DATA = ReadVergFile(DATA, strrep(DATA.layoutfile,'.layout','last.layout'));
    elseif ~isempty(layoutfile)
        DATA = ReadVergFile(DATA, layoutfile);
        fprintf('Using layout in %s\n',layoutfile);
    end
    if sum(strncmp(vpath,{'/b/bgc/matlab','/Volumes/bgc'},12)) || ~strcmp(fileparts(vpath),DATA.localmatdir)
        cprintf('red','CAREFUL!!!!!!!!!!\nVerg code is from %s not %s\n',vpath,DATA.localmatdir);
    end

    DATA = InitInterface(DATA);
    if checkforrestart
        DATA = LoadLastSettings(DATA,'interactive');
    end
    figure(DATA.toplevel);
    if ~exist('comcodes','file')
        SaveComCodes(DATA);
    end
    tt = TimeMark(tt, 'Interface');
    CheckForUpdate(DATA);
    DATA = SetExptMenus(DATA);
    SetGui(DATA);
    tt = TimeMark(tt, 'Reset Interface');
    cmdfile = ['/local/' DATA.binoc{1}.monkey '/binoccmdhistory'];
    DATA.oldcmds = scanlines(cmdfile);
    id = strncmp('Cancel',DATA.oldcmds,5);
    xid = find(strncmp('Reopen',DATA.oldcmds,5));
    if ~isempty(xid)
        lastopenstr = strrep(DATA.oldcmds{xid(end)},'Reopened ','');
    else
        lastopenstr = 'No date';
    end
    id(xid) =1;
    xid = find(strncmp('?',DATA.oldcmds,1));
    id(xid) =1;
    DATA.oldcmds = DATA.oldcmds(~id);
    DATA.oldcmds{end+1} = sprintf('History from file (%s)',lastopenstr);
    
    DATA.cmdfid = fopen(cmdfile,'a');
    myprintf(DATA.cmdfid,'Reopened %s\n',datestr(now));
    if DATA.frombinocfid > 0 || DATA.verbose(4)
        TimeMark(tt,2);
    end
    DATA = CheckStateAtStart(DATA);
    DATA.ready = 1;
    set(DATA.toplevel,'UserData',DATA);
    vpath = which('verg'); %make sure local matlab/expts is in path too
    vpath = strrep(vpath,'verg.m','expts');
    if exist(vpath,'dir')
        addpath(vpath);
    end
    if cellstrcmp('!monkeylog',vtext) || now-GetValue(DATA,'weightdate') > 30
        InterpretLine(DATA, '!monkeylog','fromverg');
    end
        
else
    DATA = get(it,'UserData');
end
end

%here parse varargin that require DATA to have been read
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{1},'close',5)
        ExitVerg(DATA);
        return;
    elseif strncmpi(varargin{j},'checkhelp',8)
        CheckCodeHelp(DATA, 'update');
    elseif strncmpi(varargin{j},'checkstim',8)
        CheckStimFile(DATA, 'update');
    elseif strncmpi(varargin{j},'setgui',5)
        ts = now;
        SetGui(DATA);
        mytoc(ts);
    elseif strncmpi(varargin{j},'getstate',5)
        ts = now;
        DATA = GetState(DATA,'commandline');
        mytoc(ts);
    elseif strcmp(varargin{j},'help')
        CodesPopup(DATA,[],'popup');
    elseif strncmpi(varargin{j},'quick',5)
        j = j+1;
        ts = now;
        stimfile = varargin{j};
        DATA.guiset = 0;
        DATA = ReadStimFile(DATA, varargin{j},'quickmenu');
        if DATA.verbose(4) mytoc(ts); end
        DATA = AddTextToGui(DATA, ['qe=' varargin{j}]);
        if DATA.verbose(4) mytoc(ts); end
        DATA = GetState(DATA,'quickexpt');
        if DATA.verbose(4) mytoc(ts); end
        SetGui(DATA,'ifneed');
        if DATA.verbose(4) mytoc(ts); end
        DATA.optionflags.afc;
        set(DATA.toplevel,'UserData',DATA);
        if DATA.verbose(4)
            fprintf('Read %s took %.3f\n',stimfile,mytoc(ts));
        end
    elseif strcmp(varargin{j},'autoreopen')
        DATA.autoreopen = 1;
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(varargin{j},'updatelog',5)
        DATA = UpdateLogFile(DATA);        
    elseif strcmp(varargin{j},'checkstart')
        DATA = CheckStateAtStart(DATA);
    elseif strncmp(varargin{j},'readstim',7)
        DATA = ReadStimFile(DATA,varargin{j+1},'-nowait');
    elseif strncmp(varargin{j},'savecom',7)
        SaveComCodes(DATA);
    elseif isfield(varargin{j},'Trials')
        DATA.Trials = varargin{j}.Trials;
        if isfield(varargin{j},'Expts')
            DATA.Expts = varargin{j}.Expts;
        end
        PlotPsych(DATA);
    end
    j = j+1;
end


if autoquit
    close(DATA.toplevel);
end


function CheckState(DATA, varargin)
% CheckState(DATA, varargin) checks values
% of params in varargin with binoc
    for j = 1:length(varargin)
        if ischar(varargin{j})
            outprintf(DATA,'%s ',varargin{j});
        end
    end
        
    if DATA.verbose(4)
        myprintf(DATA.tobinocfid,'-show','%s uf=%s bt=%.3f\n',datestr(now),DATA.binoc{1}.uf,GetValue(DATA,'bt'));
    end
    
function SaveComCodes(DATA)
    X.comcodes = DATA.comcodes;
    f = fields(DATA.helpstrs);
    for j = 1:length(X.comcodes)       
        if isfield(DATA.helpstrs,X.comcodes(j).code)
            X.comcodes(j).helpstr = DATA.helpstrs.(X.comcodes(j).code);
        end
        if isfield(DATA.helpkeys.extras,X.comcodes(j).code)
            X.comcodes(j).extra = DATA.helpkeys.extras.(X.comcodes(j).code);
        end
        if isfield(DATA.helpkeys.KeyWords,X.comcodes(j).code)
            X.comcodes(j).keys = DATA.helpkeys.KeyWords.(X.comcodes(j).code);
        end
    end
    X.optionstrings = DATA.optionstrings;
    X.helpkeys = DATA.helpkeys;
    save('comcodes','-struct','X');

        
        
    
    function ExitVerg(DATA)
    
    
    UpdateLogFile(DATA);
    SaveLayout(DATA, strrep(DATA.layoutfile,'.layout','last.layout'));

    if DATA.pipelog
        system([GetFilePath('perl') '/pipelog end']);
    end
    if isfield(DATA,'timerobj') & isvalid(DATA.timerobj)
        stop(DATA.timerobj);
    end
    for j = 2:length(DATA.windownames)
        CloseTag(DATA.windownames{j});
    end
    CloseTag(DATA.windownames{1});
    if isfigure(DATA.servofig);
        if isfield(DATA,'servotimer') & isvalid(DATA.servotimer)
            stop(DATA.servotimer);
        end
        ServoDrive('close');
    end

    fids = fopen('all');
    for j = 1:length(fids)
        name = fopen(fids(j));
        if cellstrfind(name, {'frombinoc.txt' 'binoccmdhistory' 'tobinoc'})
            if DATA.verbose(4)
                fprintf('Closing %s\n',name);
            end
            fclose(fids(j));
        else
            fprintf('Not closing file %s\n',name);
        end
    end
    WriteText(DATA.Statuslines,'/local/vergstatus.txt','backup');
    CloseChildren(DATA.toplevel);

function DATA = OpenBinocLog(DATA, type)

    if nargin == 1
        type = 'frombinoc';
    end
    
if sum(strcmp(type,{'frombinoc','both'}))    
    if DATA.frombinocfid <= 0
        DATA.frombinocfid = fopen('/local/frombinoc.txt','a');
        myprintf(DATA.frombinocfid,'Log openened %s\n',datestr(now));
    else
        fprintf('FromBinoc log already open\n');
    end
end

if sum(strcmp(type,{'tobinoc','both'}))
    if DATA.tobinocfid <= 0
        DATA.tobinocfid = fopen('/local/tobinoc.txt','a');
        myprintf(DATA.tobinocfid,'Log openened %s\n',datestr(now));
    else
        fprintf('ToBinoc log already open\n');
    end
end

function [ok, reason] = CheckDayEnd(DATA)
    ok = 0;
    reason = 'OK';
    S = DATA.binoc{1};
    if ~isfield(DATA,'Expts')
        reason = 'No Expts';
    elseif ~isfield(S,'pe') || S.Pn(1) <= 0 
        reason = 'Penetration not set';
    elseif DATA.optionflags.py
        reason = 'Human Psychophyscis';
    else
        penlog = sprintf('/local/%s/pen%d.log',S.monkey,S.Pn(1));
        nsave = 0;
        for k = 1:length(DATA.Expts)
            if isfield(DATA.Expts{k},'stored') && DATA.Expts{k}.stored > 0
                nsave = nsave+1;
            end
        end
        d = dir(penlog);
        if length(d) ==1  && now-d.datenum < 0.5 && nsave > 0 %pen log made today
            ok = 1;
        else
            reason = 'No saved Expts or penlog is old';
        end
    end
    
  


function DATA = CheckStateAtStart(DATA)
    if strcmp('NotSet',DATA.binoc{1}.ereset)
        str = 'You Can define A "reset" stimfile that is Run before loading each new Expt. Put ereset=path in verg.setup or your stimfile';
        msgbox(str);
    end
    if DATA.verbose(4) && 0
    for j = 1:length(DATA.comcodes)
       if ~isfield(DATA.helpstrs,DATA.comcodes(j).code) 
           if ~strcmp('xx',DATA.comcodes(j).code)
               fprintf('No help for %s\n',DATA.comcodes(j).code);
           end
       end
    end
    end
    
 function txt = InsertLine(txt, line, s)
 
     txt(line+1:end+1) = txt(line:end);
     txt{line} = s;

     



    


function CheckStimFile(DATA, type)
%remove out of date codes from a stim file
    txt = scanlines(DATA.stimfilename);
    goodlines = ones(size(txt));
    for j = 1:length(txt)
        [DATA, code, badcodes] = InterpretLine(DATA, txt{j}, 'test','fromCheck');
        if code < -1 
            goodlines(j) = 0;
        elseif ~isempty(badcodes)
            for k = 1:length(badcodes)
                id = regexp(badcodes{k},'[+-]')
                if length(id) > 1
                    badcode = badcodes{k}(1:id(2)-1);
                    regexprep(badcodes{k},'[a-Z]+[+-].*','$1');
                else
                    badcode = badcodes{k};
                end
                txt{j} = strrep(txt{j},badcode,'');
            end
        end
    end
txt = txt(find(goodlines));
if strcmp(type,'makenew')
outfile = strrep(DATA.stimfilename,'.stm','');
outfile = [outfile '.new'];
else
    outfile = DATA.stimfilename;
    BackupFile(outfile,'print');
end
fprintf('saving stim file to %s\n',outfile);
WriteText(txt, outfile);

function CheckCodeHelp(DATA, type)
    
    if nargin ==1
        type = 'list';
    end
    [scodes,sid] = sort({DATA.comcodes.code});
    helpfile = [DATA.localmatdir '/helpstrings.txt'];
    txt = scanlines(helpfile);
    for j = 1:length(DATA.comcodes)
        k = sid(j);
       if ~isfield(DATA.helpstrs,DATA.comcodes(k).code) 
           if ~strcmp('xx',DATA.comcodes(k).code) && ~bitand(DATA.comcodes(k).group,1024)
               fprintf('No help for %s\n',DATA.comcodes(k).code);
               if strcmp(type,'update')
               prevcode = DATA.comcodes(sid(j-1)).code;
               id = find(strncmp(prevcode,txt,length(prevcode)));
               if ~isempty(id)
                   txt = InsertLine(txt,id(1)+1,['#*' DATA.comcodes(k).code ' ' DATA.comcodes(k).label]);
               end
               end
           end
       end
    end
    
    if strcmp(type,'update')
        outfile = strrep(helpfile,'.txt','.new');
        WriteText(txt, outfile);
    end
    
  
function WriteText(txt, name, varargin)
    j = 1; 
    while j <= length(varargin)
        if strncmpi(varargin{j},'backup',5)
            BackupFile(name);
        end
        j = j+1;
    end
    
        fid = fopen(name,'w');
        if fid > 0
        for j = 1:length(txt)
            fprintf(fid,'%s\n',txt{j});
        end
        fclose(fid);
        end
        
    
function [DATA, codetype, badcodes, nlines] = InterpretLine(DATA, line, varargin)

    
    
setlist = 0;  %% don't update gui for every line read.
codetype = 0;
frombinoc = 0;
sendtobinoc =0;
paused = 0;
src = 'unknown';
srcchr = 'U';
badcodes = {};
nlines = 0;
setgui = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'from',4)
        src = varargin{j};
        if strcmp(src,'frombinoc')
            frombinoc = 1;
            srcchr = 'B';
        elseif strcmp(src,'frombinocT') %timer initiate read
            frombinoc = 2;
            srcchr = 'T';
        elseif strcmp(src,'fromgetstate')
            srcchr = 'G';
            frombinoc = 3;
        elseif strcmp(src,'fromgui')
            sendtobinoc = 1;
        elseif strcmp(src,'fromverg')
            srcchr = 'V';
        end
    elseif strncmpi(varargin{j},'tobinoc',4)
        sendtobinoc = 1;
        srcchr = 'V';
    elseif strncmpi(varargin{j},'paused',6)
        paused = 1;
    end
    j = j+1;
end
    
if nargin < 2 || isempty(line)
        return;
end

pstr = '';
if frombinoc ~= 2 
    paused = PauseRead(DATA, 1);
end
ts = now;
strs = textscan(line,'%s','delimiter','\n');
nlines = length(strs{1});
if DATA.verbose(1) && length(strs{1}) > 1
    fprintf('Reading %d lines\n',length(strs{1}));
end

for j = 1:length(strs{1})
 %   fprintf('%s\n', strs{1}{j});
    s = regexprep(strs{1}{j},'\s+\#.*$','');
%    s = strs{1}{j};
    eid = strfind(s,'=');
    codematches=0;
    codelen = 0;
    if ~isempty(eid)
        code = s(1:eid(1)-1);
        value = s(eid(1)+1:end);
        codelen = eid(1);
    else
        value = '';
        code = s;
    end
    if frombinoc && (DATA.verbose(7) || DATA.frombinocfid > 0)
        if ~strncmp(strs{1}{j},'SENDING00000',12)
            if DATA.verbose(7)
                fprintf('%s**\n',strs{1}{j});
            end
            myprintf(DATA.frombinocfid,'%s%c%d %s\n',datestr(now),srcchr,DATA.inexpt,strs{1}{j});
        elseif DATA.verbose(7) && DATA.frombinocfid %if both set, record every notification
            myprintf(DATA.frombinocfid,'%s%c%d %s\n',datestr(now),srcchr,DATA.inexpt,strs{1}{j});
        end
    end
    
    donestr = 0;
    
    if frombinoc == 0
        donestr = 1;
        oldsend = sendtobinoc;
        if sum(strncmp(s,{'!mat=' '!matnow' 'timerperiod' 'read=' 'pause' 'user=' '!onestim' 'pausetimeout'},5))
            sendtobinoc = 0;
            codetype = -1;
            if strncmp(s,'!mat',4) && ~isempty(value)
                DATA.matexpres = [];
                if strcmp(src,'fromstim')
                    DATA.matexpt = value;
                    if strncmp(s,'!matnow',7)
                        eval(value);
                    end
                else
                    fprintf('Calling %s from %s\n',value,src);
                    DATA.matexpres = binoceval(DATA, value);
                    DATA.matlabwasrun = 1;
                end
                SendCode(DATA, 'exp');
            elseif strncmp(s,'!onestim',8) 
                outprintf(DATA,'%s\n',s);
                pause(0.2);
                DATA = DrainBinocPipe(DATA,'waitforstim');
                return;
            elseif strncmp(s,'timerperiod',10) %verg special
                DATA.timerperiod = sscanf(value,'%f');
                sendtobinoc=0;
            elseif strncmp(s,'read',4) %read a set of instructions
                %            fid = fopen(value,'r');
            elseif strncmp(s,'user',4)
                estr = s(eid(1)+1:end);
                DATA.userstrings = {DATA.userstrings{:} estr};
                sendtobinoc=oldsend;
                codetype =-1; %dont send to binoc
            elseif strncmp(s,'pausetimeout=',10)
                DATA.pausetimeout = sscanf(value,'%d');
            elseif strncmp(s,'pause',5)
                if ~isempty(value)
                    DATA.readpause = str2num(value);
                elseif strncmp(s,'pauseforbinoc',9)
                    DATA = DrainBinocPipe(DATA);
                end
                %       pause(DATA.readpause);
                DATA.pausetime = now;
            end  %end of NotBinoc group
        elseif strncmp(s,'DATA',4) %lines Begining DATA are read as matlab instructions
            try
                eval(s);
                SetData(DATA);
            catch ME
                CheckExceptions(ME);
            end
            sendtobinoc = 0;
        elseif codelen > 8  % long lines from files
            sendtobinoc = 0;
            codetype = -1;
            if length(code) > 4 & sum(strcmp(code,DATA.windownames))
                iw = find(strcmp(code,DATA.windownames));
                DATA.winpos{iw} = sscanf(value,'%d');
            elseif strncmp(s,'mainwindow=',10) %old style
                DATA.winpos{1} = sscanf(value,'%d');
            elseif strncmp(s,'autoreopen=',10)
                DATA.autoreopen = sscanf(value,'%d');
                if DATA.autoreopen == 3
                    DATA.autoreopen=1;
                    DATA.autorestart = 1;
                else
                    DATA.autorestart = 0;
                end
            elseif strncmp(s,'autostart',9)
                    DATA.autorestart = 2;
            elseif strncmp(s,'optionwinpos=',10)
                DATA.winpos{2} = sscanf(value,'%d');
            elseif strncmp(s,'softoffwinpos=',10)
                DATA.winpos{3} = sscanf(value,'%d');
            elseif strncmp(s,'penlogwinpos=',10)
                DATA.winpos{4} = sscanf(value,'%d');
            else %unrecognized long codes from files are just kept in DATA.verg
                codetype = 0;
                if s(1) == '#'
                    sendtobinoc = 0;
                    donestr = 1;
                else
                    sendtobinoc = oldsend;
                    donestr = 0;
                    DATA.verg.(genvarname(code)) = value;
                end
            end
        elseif sum(strncmp(s, {'stepperxy' 'penwinxy' 'optionwinxy'}, 8))
            codetype = -1;
            sendtobinoc = 0;
        elseif strncmp(s, 'verbose=', 8)
            v = sscanf(s(9:end),'%d');
            DATA.verbose(1:length(v)) = v;
            sendtobinoc = 0;
            codetype = -1; %don't forward string as is to binoc
            SendCode(DATA,'verbose');
        elseif strncmp(s,'winpos=',7)
            DATA.winpos{1} = sscanf(value,'%d');
            codetype = -1;
            sendtobinoc = 0;
        elseif strncmp(s, 'slider', 6)
            codetype = -1;
            sendtobinoc = 0;
        elseif strncmp(s,'rw',2)
                val = sscanf(value,'%f');
                maxrw = GetValue(DATA,'maxrw');
                if  maxrw > 0 && val > maxrw
                    oldrw = DATA.binoc{1}.rw;
                    yn = confirm(sprintf('That Reward (%.2f) is bigger than your limit (%.2f). Proceed?',val(1),maxrw));
                    if yn == 0
                        sendtobinoc = 0;
                        codetype = -1;
                    else 
                        DATA.binoc{1}.rw = val;
                    end
                end            
        elseif s(1) == '!' %other commands
            if strcmp(s,'!monkeylog')
                if isfield(DATA,'toplevel')
                    MonkeyLogPopup(DATA.toplevel,[],'popup');
                end
            end
        else
            donestr = 0;
        end
        
    end
    
    
    if sendtobinoc
        tline = CheckLineForBinoc(strs{1}{j},DATA);
        outprintf(DATA,'%s\n',tline);
    end


    
    gotstr = 1; %default. set to 0 at end of elseifs
    if length(s) == 0 || donestr == 1
       
        
    elseif s(1) == '#' %defines stim code/label
        [a,b] = sscanf(s,'#%d %s');
%        a = a(1);
%        id = find(s == ' ');
%        DATA.comcodes(a).code = s(id(1)+1:id(2)-1);
%        DATA.comcodes(a).label = s(id(2)+1:end);
%        DATA.comcodes(a).const = a;
    elseif strncmp(s,'ACK:',4)
%        t = regexprep(s(5:end),'([^''])''','$1'''''); %relace ' with '' for matlab
%could use a different code if we want. ?1?
        s = strrep(s,char(9),'\n');%'\t' in c. ->\n. So that message is all on one line
        if strncmp(s,'ACK::',5)
            vergwarning(s(6:end),'newwin');
        elseif strncmp(s(5:end),' Error Reading Stimulus',10)
            if isappdata(DATA.toplevel,'WaitingForDlg') && getappdata(DATA.toplevel,'WaitingForDlg') ==0 && DATA.inexpt
                setappdata(DATA.toplevel,'WaitingForDlg',1);
                yn = gui.Dlg(['WARNING:' s(5:end)],DATA.toplevel,{'Continue' 'Cancel'});
                if strcmp(yn,'Cancel')
                    DATA = RunButton(DATA,'ForceCancel',1);
                end
                setappdata(DATA.toplevel,'WaitingForDlg',0);
            else
                fprintf('WARNING (waiting for dlg): %s\n',s(5:end));
            end
        else
            vergwarning(s(5:end));
        end
        DATA = AddStatusLine(DATA,s,3);
        DATA.lastmsg = s;

        
        
        
        
    elseif sum(strncmp(s,{'EXPTSTART' 'EXPTOVER' 'EXPTRUNNING' 'TESTOVER' 'CODE OVER' 'manexpvals'},8))
        if strncmp(s,'EXPTSTART',8)
            DATA.inexpt = 1;
            tic; PsychMenu(DATA); 
            tic; SetGui(DATA,'set'); 
            if DATA.verbose(4) fprintf('%s\n',s); end
        elseif strncmp(s,'CODE OVER',8)
            for j = 1:length(DATA.comcodes)
                if isempty(DATA.comcodes(j).code)
                    DATA.comcodes(j).code = 'xx';
                    DATA.comcodes(j).label = '';
                    DATA.comcodes(j).const = NaN;
                    DATA.comcodes(j).type = 'N';
                    DATA.comcodes(j).group = 512+1024;
                else
                    if ~isfield(DATA.helpstrs,code) && DATA.verbose(1)
                        fprintf('No help for %s\n',code);
                    end
                end
            end
        elseif strncmp(s,'EXPTRUNNING',8) %from GetState
            DATA.inexpt = 1;
        elseif strncmp(s,'EXPTOVER',8) %called at end or cancel
            exptcancel = 0;
            if DATA.inexpt %in case reopen pipes missed expt
                DATA.optionflags.do = 0;
                outprintf(DATA,'op=-do\n');
            end
            exptover = 1;
            if strcmp(s,'EXPTOVERSTATE') %just reporting state, not an event
                if DATA.inexpt == 0
                    exptover = 0; %put a check in here to resolve confusions
                end
            elseif strcmp(s,'EXPTOVERC') %expt cancelled
                exptcancel = 1;
            end
            
            wasexpt = DATA.inexpt;
            DATA.inexpt = 0;
            if DATA.nexpts > 0  && exptover %may be 0 here if verg is fired up after a crash
                DATA.Expts{DATA.nexpts}.End = now;
                DATA.Expts{DATA.nexpts}.last = length(DATA.Trials);
                if isfield(DATA.optionflags,'ts') && DATA.optionflags.ts
                    DATA.Expts{DATA.nexpts}.stored = 1;
                else
                    DATA.Expts{DATA.nexpts}.stored = 0;
                end
                if DATA.exptstoppedbyuser ==2 %cancelled
                    DATA.Expts{DATA.nexpts}.stored = -1;
                end
            end
            if exptover
                CheckTrialDurations(DATA,'EXPTOVER');
                if DATA.optionflags.py && DATA.optionflags.da && exptcancel == 0
                    CopyLog(DATA,'psych');
                    CopyLog(DATA,'serialpsych');
                end
            end

    %        tic; DATA = GetState(DATA); toc  %binoc sends state at end expt,
    %        before sending ExptOver
            tic; PsychMenu(DATA); 
            tic; SetGui(DATA,'set'); 
            ShowStatus(DATA);
            if ismember(frombinoc,[1 2]) && ~strncmp(s,'EXPTOVERSTATE',12)
                if ~wasexpt
                    myprintf(DATA.frombinocfid,'-show','Got EXPTOVER(%c), but not in Expt\n',srcchr);
                else
                    myprintf(DATA.frombinocfid,'-show','Expt over at %s\n',datestr(now));
                end
            end
            
%            DATA = DrainBinocPipe(DATA);

            if DATA.exptstoppedbyuser  
            %if user hist cancal/stop, dont repeat or move on to automatic next expt
                DATA.exptstoppedbyuser = 0;
                DATA.seqline = 0;
            elseif DATA.seqline > 0 && exptover
                myprintf(DATA.cmdfid,'-show','Sequence continuing from line %d: ',DATA.seqline);
%                if DATA.restartbinoc && wasexpt  %if inexpt ==0, may be anew restart
%                    DATA = RestartBinoc(DATA);
%                end
                DATA = ContinueSequence(DATA);
            elseif DATA.exptnextline > 0 && exptover
                DATA = ReadExptLines(DATA,{},'fromseq');
                DATA = RunButton(DATA,[],1);
            elseif DATA.rptexpts > 0 && wasexpt
                DATA.inexpt = 0;
                fprintf('Repeating Expt. %d to go\n',DATA.rptexpts)
                PauseRead(DATA, 1);
                if DATA.restartbinoc && wasexpt  %if inexpt ==0, may be anew restart
                    DATA = RestartBinoc(DATA);
                end
                myprintf(DATA.frombinocfid,'-show','%s Running Next of %d expts\n',datestr(now),DATA.rptexpts);
                outprintf(DATA,'//Nrpt is %d\n',DATA.rptexpts);
                DATA.rptexpts = DATA.rptexpts-1;
                it = findobj(allchild(DATA.toplevel),'flat','Tag','RptExpts');
                set(it,'string',num2str(DATA.rptexpts));
                DATA = GetState(DATA,'ExptRpt');
                DATA = uipause(now, DATA.binoc{1}.seqpause, 'Fixed Delay for repeats', DATA);
                PauseRead(DATA, 0);
                DATA = RunButton(DATA,[],1);
            end
            DATA.matlabwasrun = 0;
            if DATA.verbose(4) fprintf('%s:%d\n',s,DATA.inexpt); end
            DATA.newbinoc = 0;
            myprintf(DATA.frombinocfid,'INexpt Set to %d\n',DATA.inexpt);
        elseif strncmp(s,'TESTOVER',8)
           if isfield(DATA,'reopenstr') && ~isempty(DATA.reopenstr) && DATA.optionflags.do
               outprintf(DATA,'%s\n',DATA.reopenstr);
           end
        elseif strncmp(s,'manexpvals',10)
            sv = sscanf(s(11:end),'%f');
            DATA.Trial.sv(1:length(sv)) = sv;   
        end  %8 char codes

        
        
        
    
    elseif sum(strncmp(s,{'NewBinoc' 'confirm' 'exvals' 'fontsiz' 'fontname' 'layout' ...
            'localmatdir' 'netmatdir', 'oldelectrode' 'TOGGLE' 'rptexpts' 'STIMTYPE' ...
            'SENDING' 'SCODE=' 'status' 'stimdir' 'STIMC ' 'Unrecogn' 'emxpos' 'emypos' 'NOTIFYCLASH'},6))
        if strncmp(s,'NewBinoc',7)
            DATA = CheckForNewBinoc(DATA);
            DATA.newbinoc = 1;
            if DATA.optionflags.do %only do this when reopen pipes
                %                outprintf(DATA,'\\go\n');
            end
        elseif strncmp(s,'SENDING',7)
            codetype = -1;
            
        elseif strncmp(s,'confirm',7)
            yn = questdlg(s(8:end),'Update Check');
            if strcmp(yn,'Yes')
                fprintf('Confirm Yes\n');
                outprintf(DATA,'confirmpopup=1\n');
            end
        elseif strncmp(s,'NOTFIYCLASH',8)
            DATA = AddStatusLine(DATA,s, 3);
        elseif strncmp(s,'exvals',6)
            sv = sscanf(s(8:end),'%f');
            DATA.Trial.sv(1:length(sv)) = sv;
        elseif strncmp(s,'fontsiz',6)
            DATA.font.FontSize = str2double(value);
        elseif strncmp(s,'fontname',6)
            DATA.font.FontName = value;
        elseif strncmp(s,'layout',6)
            DATA.layoutfile = value;
        elseif strncmp(s,'localmatdir',6)
            DATA.localmatdir=value;
        elseif strncmp(s,'netmatdir',6)
            DATA.netmatdir=value;
        elseif strncmp(s,'oldelectrode=',10)
            estr = s(eid(1)+1:end);
            eid = find(strcmp(estr,DATA.electrodestrings));
            if isempty(eid)
                DATA.electrodestrings = {DATA.electrodestrings{:} deblank(estr)};
                DATA.electrodeid = length(DATA.electrodestrings);
            else
                DATA.electrodeid = eid(1);
            end
            DATA.binoc{1}.Electrode = estr;
        elseif strncmp(s,'TOGGLEEND',9)
            DATA = CheckToggleCodes(DATA);
        elseif strncmp(s,'TOGGLE',6)
            id = strfind(s,' ');
            cc = s(id(1)+1:id(2)-1);
            DATA.optioncodes.(cc).group = sscanf(s,'TOGGLE%d');
            if DATA.optioncodes.(cc).group == 4
                DATA.silentoption.(cc) = 1;
            end
            if ~isfield(DATA.optionflags,cc)
                DATA.optionflags.(cc) = 0;
            end
            DATA.optionstrings.(cc) = s(id(2)+1:end);
            DATA.togglecodesreceived = DATA.togglecodesreceived+1;
        elseif strncmp(s,'rptexpts',6)
            if ~isempty(value)
                DATA.rptexpts = sscanf(value,'%d');
            end
            if isempty(DATA.rptexpts)
                DATA.rptexpts = 0;
            end
        elseif strncmp(s,'STIMTYPE',6)
            id = strfind(s,' ');
            code = str2num(s(id(1)+1:id(2)-1))+1;
            DATA.stimulusnames{code} = s(id(2)+1:end);
        elseif strncmp(s,'SCODE',5)
            id = strfind(s,' ');
            icode = str2num(s(id(2)+1:id(3)-1))+1;
            label = s(id(3)+1:end);
            code = s(id(1)+1:id(2)-1);
            sid = find(strcmp(code,{DATA.strcodes.code}));
            if isempty(sid)
                sid = length(DATA.strcodes)+1;
            end
            DATA.strcodes(sid).label = label;
            DATA.strcodes(sid).icode = icode;
            DATA.strcodes(sid).code = code;
    elseif strncmp(s,'stimdir',6)
        DATA.binoc{1}.stimdir = value;
        SendCode(DATA,'stimdir');
    elseif strncmp(s,'status',6)
        sstr = s(8:end);
        if ~isempty(strfind(sstr,'Frames:')) && isfield(DATA.Trial,'rw')
            sstr = regexprep(sstr,' rw[0-9,\.]+',['$0,' num2str(DATA.Trial.rw)]);
            stype = 2;
        elseif strncmp(s,'statusS',7)
            sstr = s(9:end);
            stype = 2;
        else
            stype = 1;
        end
                
        if strncmp(s,'status=Can''t open Network',24)
            DATA.lastmsg = s(8:end);
            if DATA.errors(1) == 0
                vergwarning(DATA.lastmsg);
            end
            DATA.errors(1) = DATA.errors(1)+1;
            stype = 'errors';
        elseif sum(strncmp(s(8:end),{'No prefix'},7))
            stype = 'errors';            
        elseif sum(strncmp(s(8:end),{'Network Record'},10))
            DATA.lastmsg = s(8:end);
        elseif sum(strncmp(s(8:end),{'Expt Starting'},10))
            if ~isempty(DATA.Expts)
                sstr = ['Expt ' Expt2Name(DATA.Expts{end}) ' Starting at ' datestr(now)];
            end
            stype = 'expt';
        elseif sum(strncmp(s(8:end),{'Expt over'},9))
            if(DATA.optionflags.ts)
                savestr = 'Saved';
            else
                savestr = '(Not Saved)';
            end
            sstr = sprintf('%s %s',s(8:end),savestr);
            stype = 'expt';
        elseif sum(strncmp(s(8:end),{'Expt'},4))
                savestr = 'Saved';
            stype = 'expt';
        elseif sum(strncmp(s(8:end),{'uf='},3))
            DATA.smrfile = s(11:end);
            stype = 'expt';
            setgui = 1;
        end
        DATA = AddStatusLine(DATA,sstr,stype);
        if ishandle(DATA.statusitem)
            if ~isempty(DATA.lastmsg)
                set(GetFigure(DATA.statusitem),'Name',DATA.lastmsg);
            end
        end
        if DATA.inexpt == 0 && isfield(DATA,'toplevel')
            set(DATA.toplevel,'Name',s(8:end));
        end
        if strncmp(s,'status=Grid is',14)
            DATA = AddTextToGui(DATA,s(8:end));
        elseif strncmp(s,'status=TestStimulus',19)
            ns = sscanf(s,'status=TestStimulus %d');
            if mod(ns,20) == 0
                DATA = GetState(DATA,'test');
            end
        end
        fid = strfind(s,'Frames:');
        if ~isempty(fid) && ~isempty(DATA.Trials)
            nf = sscanf(s(fid(1)+8:end),'%d/%d (%f');
            DATA.Trials(length(DATA.Trials)).Nf = nf(1);
            DATA.Trials(length(DATA.Trials)).nf = nf(2);
            DATA.Trials(length(DATA.Trials)).dur = nf(3);
        end
        %        fprintf(s);
        elseif strncmp(s,'STIMC ',6)
            DATA.trialcounts = sscanf(s(7:end),'%f');
            if length(DATA.trialcounts) < 8
                DATA.trialcounts(8) = 0;
            elseif isfield(DATA.binoc{1},'Trw')
                rwdiff = DATA.trialcounts(8) - DATA.binoc{1}.Trw;
            end
            ShowStatus(DATA);
            if DATA.verbose(3)
                fprintf('%s\n',s);
            end
        elseif strncmp(s,'Unrecog',7)
            if DATA.verbose(4)
                fprintf('%s\n',s);
            end
            DATA = AddStatusLine(DATA,s,3);
            codetype = -2;
        end  %6 char codes
    elseif sum(strncmp(s,{'Expts' 'xyfsd' 'EDONE' 'Error'},5))
        if strncmp(s,'Expts1',6)
            DATA.extypes{1} = sscanf(s(8:end),'%d');
            DATA.extypes{1} = DATA.extypes{1}+1;
            DATA = SetExptMenus(DATA);
        elseif strncmp(s,'Expts2',6)
            DATA.extypes{2} = sscanf(s(8:end),'%d');
            DATA.extypes{2} = DATA.extypes{2}+1;
            DATA = SetExptMenus(DATA);
        elseif strncmp(s,'Expts3',6)
            DATA.extypes{3} = sscanf(s(8:end),'%d');
            DATA.extypes{3} = DATA.extypes{3}+1;
            DATA = SetExptMenus(DATA);
        elseif strncmp(s, 'EDONE', 5) %finished listing expt stims
            if isfield(DATA,'toplevel')
                it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt3StimList','style','edit');
                if length(it) == 1
                    set(it,'string',DATA.exptstimlist{3});
                end
                it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt2StimList','style','edit');
                if length(it) == 1
                    ival = get(it,'value');
                    ival = min([size(DATA.exptstimlist{2},2) ival]);
                    if ival <1
                        ival = 1;
                    end
                    set(it,'string',DATA.exptstimlist{2},'value',ival);
                end
                
                it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt1StimList','style','edit');
                if length(it) == 1
                    set(it,'string',DATA.exptstimlist{1});
                end
                B = DATA.binoc{1};
%                fprintf('EDONE at %s N %d,%d,%d\n',datestr(now),B.nt,B.n2,B.n3);
                SetExptItems(DATA);
            end
            if isfield(DATA,'exptstimlist')
                for j = 1:length(DATA.exptstimlist)
                    S = DATA.exptstimlist{j};
                    for m = 1:length(S)
                        [val,n] = sscanf(S{m},'%f');
                        if n == 1
                            DATA.expvals{j}(m) = val;
                        end
                    end
                    
                end
            end
        elseif strncmp(s,'ErrorStartExpt',10)
            DATA = AddStatusLine(DATA,s,'errors');
            if strcmp(s,'ErrorStartExpt')
                         DATA.inexpt = 0;
            end
            
        elseif strncmp(s,'Error',5)
            DATA = AddStatusLine(DATA,s,'errors');
        elseif strncmp(s,'xyfsd',5)
            x = sscanf(value,'%f');
            DATA.binoc{1}.xyfsd = x(1);
            if length(x) == 4
                DATA.showxy = x(2:4);
            end
        end %end of 5 char codes
    elseif sum(strncmp(s,{'CODE' 'cwd=' 'TRES' 'over' 'Not '},4))
        if strncmp(s,'CODE',4)
            id = strfind(s,' ');
            code = str2num(s(id(2)+1:id(3)-1))+1;
            if ismember(code,[DATA.comcodes.code])
                x = 1;
            end
            str = s(id(1)+1:id(2)-1);
            if ~strcmp(str,'xx') && sum(strcmp(str,{DATA.comcodes.code}))
                x = strcmp(str,{DATA.comcodes.code});
            end
            
            DATA.comcodes(code).label = s(id(3)+1:id(end)-2);
            DATA.comcodes(code).code = s(id(1)+1:id(2)-1);
            DATA.comcodes(code).const = code;
            DATA.comcodes(code).type = s(id(end)-1);
            is = sscanf(s(id(end)+1:end),'%d.%d');
            DATA.comcodes(code).group = is(1);
            if length(is) > 1
                DATA.comcodes(code).save = is(2);
            else
                DATA.comcodes(code).save = NaN;                
            end
            try
                DATA.codeids.(DATA.comcodes(code).code) = code; %index of codes
            catch
                fprintf('Invalid field name %s\n',str);
                DATA.codeids.xx = code; %index of codes
                DATA.comcodes(code).code = 'xx';
            end
        elseif strncmp(s,'cwd=',4)
            DATA.cwd = value;
%            DATA.binoc{1}.cwd = value;
        elseif strncmp(s,'over',4)
            DATA.over = 1;
        elseif strncmp(s,'Not ',4)
            DATA.lastline = s;
        elseif strncmp(s,'TRES',4)
            if s(6) == 'G' || s(6) == 'W'
                a = sscanf(s(7:end),'%f');
            else
                a = 0;
            end
            if(s(6) == 'G' || s(6) == 'F')
                DATA.Trial.good = 1;
            elseif(s(6) == 'W')
                DATA.Trial.good = -1;
            else
                DATA.Trial.good = 0;
            end
            if ~isnan(GetValue(DATA,'psyv'))
                DATA.Trial.psyv = GetValue(DATA,'psyv');
            end
            DATA.Trial.saved = DATA.optionflags.ts;
            DATA.Trial.RespDir = a(1);
            if length(a) > 1
                DATA.Trial.tr = a(2);
            end
            DATA.Trial.rw = DATA.currentrw;
            DATA = SetTrial(DATA, DATA.Trial);
            DATA.nt = DATA.nt+1;
            DATA.Trial.Trial = DATA.nt;
            id = findstr(s,' ');
            if length(id) > 1
                DATA.Trial.id = sscanf(s(id(2)+1:end),'%d');
            end
            if DATA.verbose(3)
                fprintf('%s t%dR%d Ex%s\n',s,DATA.nt,DATA.Trial.RespDir,sprintf(' %.3f',DATA.Trial.sv));
            end
            if isfield(DATA.Trial,'RespDir')
                %Dont plot if there are more results still in the pipeline
                if length(strs{1}) < 20 || max(find(strncmp('TRES',strs{1},4))) <= j
                    pts = now;
                    DATA = PlotPsych(DATA);
                    pstr = sprintf('Psych took %.2f',mytoc(pts));
                    
                end
            end
            
            
        end %end of 4 char codes
    elseif strncmpi(s,'exp',3)
        if strncmp(s,'exps',4)
            ex = 1;
            DATA.expts{ex} = [];
            pid = strfind(s,'+');
            for k = 1:length(pid)
                if k == length(pid)
                    x = s(pid(k)+1:end);
                else
                    x = s(pid(k)+1:pid(k+1)-1);
                end
                if length(x) > 1
                eid = find(strcmp(x,{DATA.comcodes.code}));
                DATA.expts{ex} = [DATA.expts{ex} eid];
                end
            end
            ex = 1;

        elseif strncmp(s,'expt',4)
            DATA = ReadStimFile(DATA, value, 'inread');
        elseif strncmp(s,'exp=',4)
            DATA.binoc{1}.exp = value;
        else
            gotstr = 0; %didn't process it yet
        end
    elseif sum(strncmp(s,{ 'mo=' 'pf=' 'qe=' 'nt='},3))

    if strncmp(s,'mo=fore',7)
        DATA.currentstim = 1;
        DATA.binoc{1}.mo = 'fore';
    elseif strncmp(s,'mo=back',7)
        DATA.currentstim = 2;
        DATA.binoc{1}.mo = 'back';
    elseif strncmp(s,'mo=ChoiceU',10)
        DATA.currentstim = 3;
        DATA.binoc{1}.mo = 'ChoiceU';
    elseif strncmp(s,'mo=ChoiceD',10)
        DATA.currentstim = 4;
        DATA.binoc{1}.mo = 'ChoiceD';
    elseif strncmp(s,'qe=',3)
        
        id = strfind(s,'"');
        if length(id) > 1
            submenu = s(id(1)+1:id(2)-1);
            s = s(id(end)+1:end);
        else
            submenu = '';
            s = s(4:end);
        end
        if ~strcmp(s,'NotSet')
        [a,b,c] = fileparts(s);
        b = [b c];
        if s(1) ~= '/'
            if isfield(DATA.binoc{1},'stimdir')
                s = [DATA.binoc{1}.stimdir '/' s];
            elseif isfield(DATA,'stimfilename')
                s = [fileparts(DATA.stimfilename) '/' s];
            end
        end
        id = [];
        xid = [];
        if isfield(DATA.quickexpts,'filename') %check we don't alreayd have this
            id = find(strcmp(s,{DATA.quickexpts.filename}));
            xid = find(strcmp(b,{DATA.quickexpts.name}));
        end
%when binoc echose back the qe list, stimdir instrutions from the stim file are
% lost. So only add menus if name is unique.  Might mean missing some
% qexpts when verg gets settig from binoc.
        if isempty(id) && (~sum(strcmp(src,{'frombinoc' 'frombinocT' 'fromgetstate'})) || isempty(xid))
            n = length(DATA.quickexpts)+1;
            DATA.quickexpts(n).name = b;
            DATA.quickexpts(n).filename = s;
            DATA.quickexpts(n).submenu = submenu;
            if strcmp(src,'fromgui')
                ReBuildQuickMenu(DATA);
            else
                DATA.newqe = 1;
            end
        end
        end
    elseif strncmp(s,'nt=',3)
        DATA.binoc{1}.(code) = str2num(value);
    elseif strncmp(s,'pf=',3)
        s = strrep(s,'+2a','+afc');
        if s(end) == ':'
            s = s(1:end-1);
        end
        s = [s '+'];
        id = regexp(s,'[+-]');
        xid = strfind(s,':');
        id = union(id,xid);
        f = fields(DATA.optionflags);
        stimf = fields(DATA.stimflags{1});
        nflag = 1;
        if strncmp(s,'pf=0',4) %% reset list, rather than add
            DATA.showflagseq = {};
            nflag = 1;
        else
            nflag = length(DATA.showflagseq)+1;
        end
        codetype = 1;
        for k= 1:length(id)-1
            newcode = s(id(k)+1:id(k+1)-1);
            code = find(strcmp(newcode,f));
            isshown = sum(strcmp(newcode,DATA.showflagseq));
            if s(id(k)) == ':'
                codetype = 2;
            elseif codetype == 2 %need to look for these in stimglags
                code = find(strcmp(newcode,stimf));
            elseif id(k+1) < id(k)+2
                fprintf('Code too short %s\n',s(id(k):end));
            elseif isempty(code) && DATA.togglecodesreceived == 0 && isvarname(code)

                DATA.showflags.(newcode) = 1;
                if ~isshown
                    DATA.showflagseq{nflag} = newcode;
                    nflag = nflag+1;
                end
                if DATA.togglecodesreceived
                fprintf('No Code for %s\n',s(id(k):end));
                end
            elseif isempty(code)
%                fprintf('No Code for %s\n,',s(id(k):end));                
            elseif s(id(k)) == '+' && id(k+1)
                DATA.showflags.(f{code}) = 1;
                if ~isshown
               DATA.showflagseq{nflag} = f{code};
               nflag = nflag+1;
                end
            else
                if ~isempty(code)
                    DATA.showflags.(f{code}) = 0;
                end
            end
        end
        if DATA.togglecodesreceived
            DATA.binoc{1}.pf = ShowFlagString(DATA);
        end
    end %end of 3 char codes
    elseif strncmp(s,'helpfile=',9)
        DATA = AddHelpFile(DATA,s);
    elseif sum(strncmp(s, {'et' 'e2' 'e3' 'n2' 'n3' 'em' 'm2' 'm3' 'op' 'ei' 'i2' 'i3' },2))
    if strncmp(s,'et',2)
        DATA.exptype{1} = sscanf(s,'et=%s');
        DATA.binoc{1}.(code) = value;
    elseif strncmp(s,'e2',2)
        DATA.exptype{2} = sscanf(s,'e2=%s');
        DATA.binoc{1}.(code) = value;
    elseif strncmp(s,'e3',2)
        DATA.exptype{3} = sscanf(s,'e3=%s');
        DATA.binoc{1}.(code) = value;
    elseif strncmp(s,'n2',2)
        DATA.nstim(2) = sscanf(s,'n2=%d');
        DATA.binoc{1}.(code) = str2num(value);
    elseif strncmp(s,'n3',2)
        DATA.nstim(3) = sscanf(s,'n3=%d');
        DATA.binoc{1}.(code) = str2num(value);
    elseif sum(strncmp(s,{'ei' 'i2' 'i3'},2)) %store incremenest in verg as strings. Then can have "lin/log"
        DATA.binoc{1}.(code) = value;
    elseif strncmp(s,'em',2)
        DATA.mean(1) = ReadVal(s,DATA);
        DATA.binoc{1}.em = DATA.mean(1);
    elseif strncmp(s,'m2',2)
        DATA.mean(2) = ReadVal(s,DATA);
        DATA.binoc{1}.m2 = DATA.mean(2);
    elseif strncmp(s,'m3',2)
        DATA.mean(3) = ReadVal(s,DATA);
        DATA.binoc{1}.m3 = DATA.mean(3);
    elseif strncmp(s,'op',2)
        f = fields(DATA.optionflags);
        if strncmp(s,'op=0',4) %everything else off
            for k = 1:length(f) DATA.optionflags.(f{k}) = 0; end
        end
        s = strrep(s,'+2a','+afc');
        s = strrep(s,'-2a','-afc');
        s = [s '+'];
        id = regexp(s,'[+-]');
        
        for k = 1:length(id)-1
            code = find(strcmp(s(id(k)+1:id(k+1)-1),DATA.badnames));
            if length(code) == 1
                code = find(strcmp(DATA.badreplacenames{code},f));
            else
                code = find(strncmp(s(id(k)+1:id(k+1)-1),f,id(k+1)-id(k)-1));
                for j = 1:length(code)
                    if length(f{code(j)}) == id(k+1)-id(k)-1
                        code = code(j);
                        break;
                    end
                end
            end
            
            if isempty(code) 
                if DATA.togglecodesreceived
                    badcodes{end+1} = s(id(k):end-1);
                    fprintf('No Code for %s\n',badcodes{end});
                else
                    myprintf(DATA.frombinocfid,'-show','Cant set Code %s - have not received list from binoc\n',s(id(k):end-1));                    
                end
            elseif s(id(k)) == '+'
            DATA.optionflags.(f{code}) = 1;
            else
            DATA.optionflags.(f{code}) = 0;
            end
        end
    end
    elseif regexp(s,'^ch[0-9]')

    if strncmp(s,'ch10',4) %Right eye XY
        DATA.showxy(1) = strfind('-+',s(5))-1;
        id = strfind(s,'fs');
        if length(id) == 1
            DATA.binoc{1}.xyfsd = sscanf(s(id(1)+2:end),'%f');
        end
    elseif strncmp(s,'ch11',4) %L eye XY
        DATA.showxy(2) = strfind('-+',s(5))-1;
    elseif strncmp(s,'ch12',4) %binoc XY
        DATA.showxy(3) = strfind('-+',s(5))-1;
    end

    elseif strncmp(s, 'st=', 3)
        id = find(strcmp(deblank(s(4:end)),DATA.stimulusnames));
        if length(id) == 1
        DATA.stimtype(DATA.currentstim) = id;
        DATA.binocstr.st = deblank(s(4:end));
        DATA.binoc{DATA.currentstim}.st = deblank(s(4:end));
        end
      %strcodes not used any more
%    elseif sum(strcmp(code,{DATA.strcodes.code}))
%        id = strfind(s,'=');
%        if id
%            sid = find(strcmp(code,{DATA.strcodes.code}));
%            if isempty(sid)
%                DATA.binocstr.(code)=s(id(1)+1:end);
%            else
%                DATA.binocstr.(DATA.strcodes(sid).code)=s(id(1)+1:end);
%           end
%       end
    elseif strncmp(s, 'Bs', 2)
             DATA.stimtype(2) = find(strcmp(s(4:end),DATA.stimulusnames));
             DATA.binoc{1}.Bs = DATA.stimulusnames{DATA.stimtype(2)};
    else
        gotstr = 0;
    end
    if gotstr == 0
        codematches = strcmp(code,{DATA.comcodes.code});
    end
        
    if gotstr ==1 %already processed
    elseif sum(codematches)
        cid = find(codematches);
        code = DATA.comcodes(cid(1)).code;
        id = strfind(s,'=');
        if frombinoc && strcmp(code,'uf')
            if DATA.state.query
                CheckState(DATA,s, 'query');
            else
                CheckState(DATA,s);
            end
        end
        if id
%some variables are numeric internally, but best recorded via the matching code        
            if sum(strcmp(code,{'hxtype'}))
                vtype = 'C';
            else
                vtype = DATA.comcodes(cid(1)).type;
            end
            if strcmp(code,'rw')
                val = sscanf(s(id(1)+1:end),'%f (%f)');
                if length(val) > 1
                    DATA.currentrw = val(2);
                end
            end
            if  vtype == 'C'
                DATA.binoc{DATA.currentstim}.(code) = s(id(1)+1:end);
            else
                val = sscanf(s(id(1)+1:end),'%f');
                DATA.binoc{DATA.currentstim}.(code) = val;
            end
            codetype = DATA.comcodes(cid(1)).group;
            DATA = SetCode(DATA,code);
        end    
    elseif regexp(s,'^E[A-C][0-9,\-,C]') %don't let this catch other codes starting with E
        if strncmp(s,'EBCLEAR',5)
            DATA.exptstimlist{2} = {};
            DATA.nextras(2) = 0;
            DATA = CheckCustomStim(DATA,s,2);
        elseif strncmp(s,'ECCLEAR',5)
            DATA.exptstimlist{3} = {};
            DATA.nextras(3) = 0;
            DATA = CheckCustomStim(DATA,s,3);
        elseif strncmp(s,'EACLEAR',5)
            DATA.exptstimlist{1} = {};
            DATA.nextras(1) = 0;
            DATA = CheckCustomStim(DATA,s,1);
        elseif s(2) == 'C'
            n = sscanf(s(3:end),'%d');
            id = findstr(s,'=');
            if length(n)
                DATA.exptstimlist{3}{n(1)+1} = s(id(1)+1:end);
                if isfield(DATA,'toplevel')  && setlist
                    it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt3StimList');
                    if length(it) == 1
                        set(it,'string',DATA.exptstimlist{3});
                    end
                end
            end
        elseif s(2) == 'B'
            n = sscanf(s(3:end),'%d');
            id = findstr(s,'=');
            if length(n)
                DATA.exptstimlist{2}{n(1)+1} = s(id(1)+1:end);
%                fprintf(s);
                if isfield(DATA,'toplevel') && setlist
                    it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt2StimList','style','edit');
                    if length(it) == 1
                        ival = get(it,'value');
                        ival = min([size(DATA.exptstimlist{2},2) ival]);
                        set(it,'string',DATA.exptstimlist{2},'value',ival);
                    end
                end
            end
        elseif s(2) == 'A'
            n = sscanf(s(3:end),'%d');
            if n < 0 && isempty(DATA.exptstimlist{1})
                DATA.nextras(1) = abs(n(1));
            end
            id = findstr(s,'=');
            if length(n)
                DATA.exptstimlist{1}{n(1)+1+DATA.nextras(1)} = s(id(1)+1:end);
                if isfield(DATA,'toplevel') && setlist
                    it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt1StimList');
                    if length(it) == 1
                        set(it,'string',DATA.exptstimlist{1});
                    end
                end
            end
            if strncmp(s,'EACLEAR',5)
                DATA.exptstimlist{1} = {};
            end
        else
            
        end
    else
        id = strfind(s,'=');
        if id & s(1) ~= '!'
% If we get here, its a command not recognized byt binoc or by lines above.
%In general, Do NOT add these to binoc{1} - it allows the user to create
%abritrary fields there accidentally. Put then in binoc{5} (not send to binoc)
%However, catch codes that might be
%set in verg.setup before binoc has been started. 
            code = s(1:id(1)-1);
            code = strrep(code, 'electrode','Electrode');
            if isfield(DATA,'matexpres') && isfield(DATA.matexpres,'expvars') && sum(strcmp(code,DATA.matexpres.expvars))
                DATA.Trial.(code) = str2num(value);
                %this is OK 
            elseif isvarname(code) %legal name
                % was  isempty(find(strcmp(code, {'1t' '2t' '3t' '4t' ''}))) %illegal names
                if sum(strcmp(code,{'ereset'}))
                    bid = DATA.currentstim;
                else
                    bid = 5;
                end
                code = deblank(code);
                val = sscanf(s(id(1)+1:end),'%f');
                if ~isempty(val)
                    DATA.binoc{bid}.(code) = val;
                else
                    DATA.binoc{bid}.(code) = s(id(1)+1:end);
                end
                if 1 && DATA.togglecodesreceived
                    codetype = -2;
                    codematches = strcmp(code,{DATA.comcodes.code});
                    xid = find(bitand([DATA.comcodes.group],2048));
                    for k = 1:length(xid)
                        if strncmp(code,DATA.comcodes(xid(k)).code,length(DATA.comcodes(xid(k)).code))
                            codetype = xid(k);
                            codematches = strncmp(code,{DATA.comcodes.code},length(DATA.comcodes(xid(k)).code));
                        end
                    end
                    if codetype < 0 && srcchr ~= 'V'
                        s = GetValue(DATA,'expvars');
                        if strfind(s,code)                            
                            if DATA.trialcounts(6) == 0
                                fprintf('%s:Seem to have manual Code %s \n',datestr(now),code);
                            end
                        else
                            fprintf('%s:Code %s not in comcodes\n',datestr(now),code);
                        end                        
                    end
                end
                SetCode(DATA,code);
            else
                cprintf('red','Cannot Interpret %s: %s\n',src,deblank(strs{1}{j}));
                codetype = -2;
            end
        else
            if strncmp(s,'cy',2) %obsolete codes
                codetype = -2;
            end
        end
    end
    if frombinoc && ~isempty(s) && codetype >= 0
            DATA.lastline = s;
    end
end

dur = mytoc(ts);
if length(strs{1}) > 20
    myprintf(DATA.frombinocfid,'Reading %d lines took %.2f%s\n',length(strs{1}),dur,pstr);
end
if dur > 1
    fprintf('Reading %d lines took %.2f%s at %s\n',length(strs{1}),dur,pstr,datestr(now));
end
if setgui
    SetGui(DATA);
end
if frombinoc ~= 2 && paused == 0 %wasnt paused before this call.
    PauseRead(DATA,0);
end

function DATA = AddStatusLine(DATA, str, type)
%add a statusline    
    %type 1 = regular status linte
    %2 = Stimulus lines
    %3 = Errors/warnings
    %4 = Expt Control
    DATA.Statuslines{end+1} = str;
    if ischar(type)
        type = find(strcmp(type,{'status' 'trial' 'errors'  'expt' 'comment'}));
        if isempty(type)
            type = 1;
        end
    end
    DATA.statustypes(length(DATA.Statuslines)) = type;
    if isfield(DATA.showstatus,'update') && DATA.showstatus.update
    ShowStatusStrings(DATA);
    end

        
function QueryBinoc(DATA,code, varargin);

function DATA = SetCode(DATA,code)
%DATA = SetCode(DATA,val)
% Does additional steps needed for some codes, like 'Electrode'

if DATA.currentstim == 1 && sum(strcmp(code,{'Electrode'}))
%don't call GetValue unles needed - wastest time
    val = GetValue(DATA,code);
    if strcmp(code,'Electrode')
        id = find(strcmp(val,DATA.electrodestrings));
        if isempty(id)
            DATA.electrodestrings{end+1} = deblank(val);
            DATA.electrodeid = length(DATA.electrodestrings);
        else
            DATA.electrodeid = id;
        end
    end
end


function DATA = CheckCustomStim(DATA, s, n)
    if s(end) == '*'
        DATA.customstimlist(n) = 1;
    else
        DATA.customstimlist(n) = 0;
    end
    
function s = AddCustomStim(DATA, s, n)
    
    pre = {'EA' 'EB' 'EC'};
    for k = 1:length(n)
    o = 1+DATA.nextras(n(k));
    
    if DATA.customstimlist(n(k))
        for j = o:length(DATA.exptstimlist{n(k)})
            s = [s sprintf('%s%d=%s\n',pre{n(k)},j-o,DATA.exptstimlist{n(k)}{j})];
        end
    end
    end
    

function [code, codeid] = FindCode(DATA, s)
    code = '';
    codeid = 0;
    id = strfind(s,'=');
    if id
        code = s(1:id(1)-1);
        codeid = find(strcmp(code,{DATA.comcodes.code}));
    elseif length(s) > 1
        id = find(strcmp(s(1:2),{DATA.comcodes.code}));
        codeid = id;
        if length(id) == 1
            code = DATA.comcodes(id).code;
        end
    end
        
function str =  ShowFlagString(DATA)    
    str = 'pf=';
    for j = 1:length(DATA.showflagseq)
        f = DATA.showflagseq{j};
        if ~isfield(DATA.showflags,f)
            fprintf('Showflags no field %s\n',f);
        elseif DATA.showflags.(f)
            str = [str '+' f];
        end
    end
    
    
function res = binoceval(DATA, str)
    id = strfind(str,'$');
    while ~isempty(id)
        sid = regexp(str(id(1)+1:end),'[\s\,\)]')+id(1);
        f = str(id(1)+1:sid(1)-1);
        if isfield(DATA.matexpvars,f)
            str = [str(1:id(1)-1) num2str(DATA.matexpvars.(f)) str(sid(1):end)];
        elseif isfield(DATA.binoc{1},f)
            str = [str(1:id(1)-1) num2str(DATA.binoc{1}.(f)) str(sid(1):end)];
        end
        id = strfind(str,'$');
    end
    h = waitbar(0.5,sprintf('Running %s',strrep(str,'_',' ')));
    fprintf('Running %s\n',str);
    outprintf(DATA,'#!mat=%s %s\n',datestr(now),str);
    res = eval(str);
    if isfield(res,'version')
        outprintf(DATA,'stimver=%.3f\n',res.version);
    end
    delete(h);
    
    
function DATA = CheckToggleCodes(DATA)    
    f = fields(DATA.showflags);
    of = fields(DATA.optionflags);
    
    for j = 1:length(f)
        if sum(strcmp(f{j},of)) == 0
            fprintf('Option %s Not supported\n',f{j});
            DATA.showflags = rmfield(DATA.showflags,f{j});
            id = find(strcmp(f{j},DATA.showflagseq));
            id = setdiff(1:length(DATA.showflagseq),id);
            DATA.showflagseq = DATA.showflagseq(id);
        end
    end

function [DATA, details] = ReadExptLines(DATA, strs, src,varargin)

    if isempty(strs)
        strs = DATA.explines
        firstline = 1;
    else
        firstline = 1+DATA.exptnextline;
    end
    badcodes = 0;
    [a,b] = cellstrcmp('-nowait',varargin);
    if a
        oargs = varargin(find(b));
    else
        oargs = {};
    end
    
    for j = firstline:length(strs)
        tline = strs{j};
        if strncmp(tline,'next',4)
            DATA.exptnextline = j;
            break;
        end
        if ~ischar(DATA.binoc{1}.monkey)
            fprintf('Monkeyame corrupted (%s)\n',tline);
        elseif ~isempty(strfind(tline,'$MNK')) && ~isempty(DATA.binoc{1}.monkey)
            tline = strrep(tline,'$MNK',DATA.binoc{1}.monkey);
        end
        monkey = DATA.binoc{1}.monkey;
        if strcmp(tline,'immode=preload')
            fprintf('Substituting imload for immode\n');
            tline = 'imload=preload';
        end
        if strncmp(tline,'openbinoclog',12)
            DATA = OpenBinocLog(DATA,'frombinoc');
            if strncmp(tline,'openbinoclogs',13)
                DATA = OpenBinocLog(DATA,'tobinoc');
            end
        else
            [DATA, type] = InterpretLine(DATA,tline, src);
        end
        if DATA.perfmonitor
            myprintf(DATA.frombinocfid,'%.3f file %s\n',mytoc(DATA.starttime),tline);
        end
  %if filename set before monkey, set monkeyname first
        if strncmp(tline,'uf=',3) && strcmp(DATA.binoc{1}.monkey,'none')
            monkey = GetMonkeyName(DATA.binoc{1}.uf);
            if ~isempty(monkey)
                DATA.binoc{1}.monkey = monkey;
                SendCode(DATA,'monkey');
            else
                myprintf(DATA.frombinocfid,'Can''t find Monkey Name in %s\n',tline)
            end
        end
        if strncmp(tline,'monkey',6)
            if ~strcmp(monkey,DATA.binoc{1}.monkey) %changed Subject
                NewMonkey(DATA, monkey);
            end
        end
        if DATA.over
            DATA.overcmds = {DATA.overcmds{:} tline};
        else
            if type >= 0 %don't send lines type < 0to binoc 
                tline = CheckLineForBinoc(tline, DATA);
%                outprintf(DATA,'%s\n',tline,'-nowait');
                outprintf(DATA,'%s\n',tline,oargs{:});
            else
                if DATA.verbose(4)
                    fprintf('%s Not sent\n',tline);
                end                 
            end
        end
        if type == -2
            badcodes = badcodes+1;
        end
    end
    if j >= length(DATA.exptlines)
            DATA.exptnextline = 0;
    end
    outprintf(DATA,'\neventcontinue\nEDONE\n');
    DATA = ReadFromBinoc(DATA,'expect');
    for ex = 1:3
        if length(DATA.expts{ex})
            id = find(~ismember(DATA.expmenuvals{ex}, DATA.expts{ex}));
            DATA.expmenuvals{ex} = [DATA.expts{ex} DATA.expmenuvals{ex}(id)];
        end
        DATA = SetExptMenus(DATA);
    end
    details.badcodes = badcodes;

function vergwarning(s, varargin)

    persistent lastmsg;
    persistent lastcalltime;
    persistent showpopup;
    newwindow = 0;
toconsole = 1;
tellbinoc = 0;
    if isempty(showpopup)
        showpopup = 1;
    end
    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'newwin',3)
            newwindow = 1;
        elseif strncmpi(varargin{j},'tellbinoc',6)
            tellbinoc = 1;
        end
        j = j+1;
    end
    if strcmp(s,'-nopopup')
        showpopup = 0;
        return;
    elseif strcmp(s,'-popup')
        showpopup = 1;
        return;
    end
   OldPos = [];
   ts = now;
   
   if strncmp(s,lastmsg,10) %repeated message
       secdiff = (ts - lastcalltime) * (24 * 60 * 60);
       if secdiff < 1
           lastcalltime = ts;
           return;
       elseif strncmp(s,' Trial 0 Initial Serial Line Not Responding',30) && secdiff < 20
           lastcalltime = ts;
           return;
       end
   end  
   beep;
   CreateStruct.Interpreter = 'tex';
   if newwindow == 0
        OldFig = findobj(allchild(0),'flat','tag','Msgbox_Binoc Warning','Name','Binoc Warning');
        if ~isempty(OldFig)
           OldPos = get(OldFig(1),'position');
        end
        CreateStruct.WindowStyle='replace';
   else
        CreateStruct.WindowStyle='non-modal';
   end
   try
       
       if tellbinoc
           outprintf(DATA,'#%s\n',s)
       end
       if showpopup
           h = msgbox(split(s,'\\n'),'Binoc Warning','warn',CreateStruct);
           p = ScaleWindow(h,2);
           if ~isempty(OldPos)
               set(h,'Position', [OldPos(1:2) p(3:4)]);
               drawnow;
           end
       end
       if toconsole
        fprintf('WARNING: %s at %s\n',s,datestr(now));
       end
   end
   
   lastmsg = s;
   lastcalltime = ts;
   
function line = CheckLineForBinoc(tline, DATA)
    if strncmp(tline,'op',2)
        tline = strrep(tline,'+2a','+afc');
    end
% more recent Spike2 seems not to need this
%But matlab needs it for writing to serial line
    if strncmp(tline,'uf',2)
%        tline = regexprep(tline,'\\\','\'); %dont go '\\' -> '\\\\'
%        tline = regexprep(tline,'\\','\'); %dont go '\\' -> '\\\\'
%        tline = regexprep(tline,'\','\\');
%beware fprintf('%s',tline) and fprintf(tline) behave differetnly with \\
        tline = regexprep(tline,'/','\\'); %Spike 2 needs \ not /
        if ~isempty(DATA.pcdrive)
            tline = strrep(tline,'$D',DATA.pcdrive);
        end
    else
        tline = regexprep(tline,'\','/');
    end
    tline = regexprep(tline,'\s+\#.*\n','\n'); %remove comments
    line = strrep(tline,'\s+\n','\n');
    
    
    
function DATA = ReadSetupFile(DATA, name, varargin)
fid = fopen(name,'r');
if fid > 0
    tline = fgets(fid);
    while ischar(tline)
        id = strfind(tline,'=');
        value = '';
        if ~isempty(id)
            value = tline(id+1:end);
        end
        if strncmpi(tline,'electrode=',8) && ~isempty(value)
            DATA.electrodestrings{1+length(DATA.electrodestrings)} = deblank(value);
        end
        if strncmp(tline,'tty2=',5) && ~isempty(value)
            DATA.servoport = deblank(value);
        end
        tline = fgets(fid);
    end
    fclose(fid);
end

function [DATA, details] = ReadStimFile(DATA, name, varargin)
   
    setall = 0;
    inread = 0;
    src = 'unknown';
    details.badcodes = 0;
    linesrc = 'fromstim';
    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'inread',4)
            inread = 1;
        elseif strncmpi(varargin{j},'init',4)
            setall = 1;
        elseif strncmpi(varargin{j},'fromverg',7)
            linesrc = 'fromverg';
        elseif strncmpi(varargin{j},'quickmenu',8)
            src = varargin{j};
        end
        j = j+1;
    end
%if this is the first load since running an expt, Call the reset
%file
DATA.newqe = 0;
if DATA.newexptdef == 0 && isfield(DATA.binoc{1},'ereset') && ~strcmp('NotSet',DATA.binoc{1}.ereset);
    myprintf(DATA.cmdfid,'-show','#Resetting Expt with %s\n',DATA.binoc{1}.ereset);
    DATA.newexptdef = 1;
    outprintf(DATA,'newexpt\n');
    DATA = ReadStimFile(DATA, DATA.binoc{1}.ereset);
end
outprintf(DATA,'//qe%s\n',name);

fid = fopen(name,'r');
if fid > 0
    if strcmp(src,'quickmenu') %undo some things whenever load new quick
        DATA.matexpt = []; 
    end
        
            
    if  inread == 0
        outprintf(DATA,'\neventpause\n');
    end
    for j = 2:length(DATA.overcmds) %commands to execute at and
            outprintf(DATA,[DATA.overcmds{j} '\n']);
        end
    DATA.over = 0;
    DATA.overcmds = {};
    DATA.exptnextline = 0;
    a = textscan(fid,'%s','delimiter','\n');
    DATA.exptlines = a{1};
    fclose(fid);
    sid = find(strncmp('sequence',DATA.exptlines,8));
    if ~isempty(sid)
        popup = 'popup';
        seqlines = DATA.exptlines(sid+1:end);
        if sid > 1
            if strcmp(DATA.exptlines{sid},'sequencenew')
                popup = 'popupnew';
            end
%if sequcence appears after line 1, only send lines about this to binoc
            DATA.exptlines = DATA.exptlines(1:sid-1);
        end
    end
    if isempty(sid) || sid > 1
        [DATA, details] = ReadExptLines(DATA,DATA.exptlines,linesrc,varargin{:});
        if details.badcodes > 1
            fprintf('Choose Fix from the file menu, or run verg([],''checkstim'') to remove bad/old codes from %s\n',name);
        end
    else %just a sequence file
        outprintf(DATA,'\neventcontinue\n');
    end
    if ~isempty(sid)
        SequencePopup(DATA,seqlines,popup);
    end
else
    try
    msgbox(sprintf('Can''t read stimfile %s',name),'Stimfile Read Error','error');
    catch
        fprintf('Can''t read stimfile %s\n',name)
    end
end
if strcmp(src,'quickmenu') && DATA.newqe > 0
    ReBuildQuickMenu(DATA);
end

if setall
    DATA = ReadVergFile(DATA, DATA.layoutfile);
end
DATA.newexptdef = 1;
DATA.userstrings = unique(DATA.userstrings);

function DATA = ReadVergFile(DATA, name, varargin)
%Read file htat only affects verg, not binoc (e.g. layout)
%so these lines are not send on to binoc   
    setall = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'init',4)
            setall = 1;
        end
        j = j+1;
    end
        
fid = fopen(name,'r');
if fid > 0
    tline = fgets(fid);
    while ischar(tline)
        DATA = InterpretLine(DATA,tline,'vergfile');
        tline = fgets(fid);
    end
    fclose(fid);
else
    msgbox(sprintf('Can''t read Verg file %s',name),'VergFile Read Error','error');
end




function SendState(DATA, varargin)

    sendview = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'All',4)
            sendview = 1;
        end
        j = j+1;
    end
    
    f = fields(DATA.binoc{1});
    
    if ~isfield(DATA,'codeids') %not connected to binoc
        return;
    end
    outprintf(DATA,'\neventpause\n');
    
    
    outprintf(DATA,'mo=fore\n');
    if sendview
        SendCode(DATA,{'px' 'py' 'vd'});
    end
    SendCode(DATA,{'monkey'});
    f = fields(DATA.binocstr);
    for j = 1:length(f)
        if length(DATA.binocstr.(f{j})) > 0
            outprintf(DATA,'%s=%s\n',f{j},DATA.binocstr.(f{j}));
        end
    end
    outprintf(DATA,'uf=%s\n',DATA.binoc{1}.uf);
    
    if DATA.showxy(3)
        outprintf(DATA,'ch12+\n'); %R
    end
    if DATA.showxy(2)
        outprintf(DATA,'ch11+\n'); %L
    end
    if DATA.showxy(1)
        outprintf(DATA,'ch10+\n'); %B
    end
%send op first, as it may change the way what follows is interpreted
% eg if op=+bf (back fixed)
    SendCode(DATA,{'optionflag'}); 

    if length(DATA.binoc) > 1 && isstruct(DATA.binoc{2})
    f = fields(DATA.binoc{2});
    outprintf(DATA,'mo=back\n');
    outprintf(DATA,'st=%s\n',DATA.stimulusnames{DATA.stimtype(2)});

    for j = 1:length(f)
        if sum(strcmp(f{j},DATA.redundantcodes))
        elseif ischar(DATA.binoc{2}.(f{j}))
        outprintf(DATA,'%s=%s\n',f{j},DATA.binoc{2}.(f{j}));
        else
        outprintf(DATA,'%s=%.6f\n',f{j},DATA.binoc{2}.(f{j}));
        end
    end
    end
    SendChoiceTargets(0,DATA);
    outprintf(DATA,'mo=fore\n');
    outprintf(DATA,'st=%s\n',DATA.stimulusnames{DATA.stimtype(1)});
    f = fields(DATA.binoc{1});
    for j = 1:length(f)
        if isfield(DATA.codeids,f{j})
            cid = DATA.codeids.(f{j});
        else
            cid = 1; 
        end
        if DATA.comcodes(cid).group > 511
            if DATA.verbose(4)
            fprintf('Not sending %s\n',f{j});
            end
        elseif sum(strcmp(f{j},DATA.redundantcodes))
            if DATA.verbose(4)
                fprintf('Not sending %s\n',f{j});
            end
        elseif isempty(DATA.binoc{1}.(f{j}))
            if DATA.verbose(4)
                fprintf('Not sending Empty value %s\n',f{j});
            end
        elseif ischar(DATA.binoc{1}.(f{j}))
            outprintf(DATA,'%s=%s\n',f{j},DATA.binoc{1}.(f{j}));
        else
            outprintf(DATA,'%s=%s\n',f{j},sprintf('%.6f ',DATA.binoc{1}.(f{j})));
        end
    end
    outprintf(DATA,'clearquick\n');
    for j = 1:length(DATA.quickexpts)
        if length(DATA.quickexpts(j).submenu)
            outprintf(DATA,'qe="%s"%s\n',DATA.quickexpts(j).submenu,DATA.quickexpts(j).filename);
        else
            outprintf(DATA,'qe=%s\n',DATA.quickexpts(j).filename);
        end
    end
%some special codes are sent to binoc differently.  
%Codes that go wrong during ReopenPipes probably need to be added here.
    SendCode(DATA,'vve');
    SendCode(DATA,'optionflag');
    SendCode(DATA,'expts');
    SendCode(DATA,'ed'); %need to tell binoc !seted=
    if DATA.electrodeid > 0
        outprintf(DATA,'Electrode=%s\n',DATA.electrodestrings{DATA.electrodeid});
    end
        
    outprintf(DATA,'\neventcontinue\n');

function WriteStr(DATA, fid, varargin)
    str = sprintf(varargin{:});
    if fid <= 0
        outprintf(DATA, str);
    else
        fprintf(fid, str);
    end
    
    
    
function SendChoiceTargets(fid, DATA)
    if length(DATA.stimtype) < 3
        return;
    end
    if DATA.stimtype(3)
        WriteStr(DATA, fid,'mo=ChoiceU\n');
        WriteStr(DATA, fid,'st=%s\n',DATA.stimulusnames{DATA.stimtype(3)});
        f = fields(DATA.binoc{3});
        for j = 1:length(f)
            [s, lbl] = CodeText(DATA, f{j},'ChoiceU');
            WriteStr(DATA, fid,'%s\t#ChoiceU %s\n',s,lbl);
        end
%        StimToggleString(DATA,3);
    end
    if DATA.stimtype(4)
        WriteStr(DATA, fid,'mo=ChoiceD\n');
        WriteStr(DATA, fid,'st=%s\n',DATA.stimulusnames{DATA.stimtype(4)});
        f = fields(DATA.binoc{4});
        for j = 1:length(f)
            [s, lbl] = CodeText(DATA, f{j},'ChoiceD');
            WriteStr(DATA, fid,'%s\t#ChoiceD %s\n',s,lbl);
        end
 %       StimToggleString(DATA,3);
    end
    

function SaveExpt(DATA, name)
    bname = name;
    if strfind(name, '.stm')
        bname = strrep(name,'.stm','.bstm');
    else
        bname = [name '.bstm'];
    end
    BackupFile(name,'print');
    outprintf(DATA,'\\savefile=%s\n',bname);
    fid = fopen(name,'w');

    fprintf(fid,'mo=fore\n');
    fprintf(fid,'st=%s\n',DATA.stimulusnames{DATA.stimtype(1)});
    fprintf(fid,'%s\n',CodeText(DATA, 'expts'));
    fprintf(fid,'%s\n',CodeText(DATA, 'nr'));
    fprintf(fid,'uf=%s\n',DATA.binoc{1}.uf);
    
    fprintf(fid,'mo=back\n');
    fprintf(fid,'st=%s\n',DATA.stimulusnames{DATA.stimtype(2)});
    f = fields(DATA.binoc{2});
    for j = 1:length(f)
        [s, lbl] = CodeText(DATA, f{j},'back');
        fprintf(fid,'%s\t#Back %s\n',s,lbl);
    end
    SendChoiceTargets(fid, DATA);
    fprintf(fid,'mo=fore\n');
    fprintf(fid,'st=%s\n',DATA.stimulusnames{DATA.stimtype(1)});
%order of code wrting set by order of fields
    allf = fields(DATA.binoc{1});
    [a,b] = ismember(allf,{DATA.comcodes.code});
    f = {DATA.comcodes(sort(b)).code};
    for j = 1:length(f)
        [s, lbl, type, cid] = CodeText(DATA, f{j});
        if bitand(type, 512) == 0 && DATA.comcodes(cid).save == 1;
            fprintf(fid,'%s\t#%s\n',s,lbl);
        elseif strcmp(GetValue(DATA,'ui'),'bgc')
            fprintf('%s not saved %s\n',s,lbl);
        end        
    end
    f = fields(DATA.binocstr);
    for j = 1:length(f)
        id = find(strcmp(f{j},{DATA.comcodes.code}));
        if length(id) == 1
            fprintf(fid,'%s=%s\t#%s\n', f{j}, DATA.binocstr.(f{j}),DATA.comcodes(id).label);
        else
            fprintf(fid,'#Error writing %s\n',f{j});
        end
    end
    fprintf(fid,'%s\n',CodeText(DATA, 'optionflag'));
    fprintf(fid,'%s\n',CodeText(DATA, 'pf'));
    for j = 1:length(DATA.quickexpts)
        if isempty(DATA.quickexpts(j).submenu)
        fprintf(fid,'qe=%s\n',DATA.quickexpts(j).filename);
        else
        fprintf(fid,'qe="%s"%s\n',DATA.quickexpts(j).submenu,DATA.quickexpts(j).filename);
        end
    end
    for j = 1:length(DATA.helpfiles)
        fprintf(fid,'helpfile="%s"%s\n',DATA.helpfiles(j).label,DATA.helpfiles(j).filename);
    end
    fclose(fid);
    fprintf('Saved %s\n',name);

function DATA = OpenPipes(DATA, readflag)
        
DATA.outpipe = '/tmp/binocinputpipe';
DATA.inpipe = '/tmp/binocoutputpipe';

rbusy = 0;

[a, pstr] = system('ps -e | grep binoclean');
if isempty(strfind(pstr, 'binoclean.app'))
    if DATA.autorestart
        if DATA.autorestart == 2 %only start at beginning
            DATA.autorestart = 0;
        end
        [a,b] = system('open /local/bin/binoclean.app');        
        pause(1);
        binocisrunning = ~a; 
        if a ~= 0
            vergwarning('Tried to Open Binoc But Failed!!!!\n');
        end
    else
        fprintf('Binoc is Not Running');
        binocisrunning = 0;
    end
else
    binocisrunning = 1;
end

if DATA.network
    DATA.binocisup = 1;
    if binocisrunning == 0 && DATA.autorestart
        fprintf('Please Wait while I start binoc....\n');
        [DATA.binoncisup,b] = system('open /local/bin/binoclean.app');        
        if DATA.binocisup %launch succeeded. Wait until its ready
            d = [];
            while isempty(d)
             d = dir('/tmp/binocisnew');
             pause(0.1);
            end
            fprintf('Thank you. Binoc is running now\n');       
        end
    end
    warning('off','MATLAB:urlread:ReplacingSpaces');
else
    
    if DATA.outid > 0
        fclose(DATA.outid);
    end
    if DATA.inid > 0
        fclose(DATA.inid);
    end
    
    if exist(DATA.outpipe)  && binocisrunning
        %if binoc crashed out leaving pipes behind, this will freeze.
        DATA.outid = fopen(DATA.outpipe,'w');
    else
        DATA.outid = 1;
    end
    if exist(DATA.inpipe) && binocisrunning
        DATA.inid = fopen(DATA.inpipe,'r');
    else
        DATA.inid = 1;
    end
end
outprintf(DATA,'NewMatlab\n');
ts = now;
[DATA, str] = ReadFromBinoc(DATA,'reset');
n = 1;
while DATA.togglecodesreceived == 0 && mytoc(ts) < 4
    pause(0.01);
    if binocisrunning == 0
        outprintf(DATA,'NewMatlab\n');
    end
    [DATA, str] = ReadFromBinoc(DATA);
    n = n+1;
end
DATA.readdur = mytoc(ts);
if DATA.verbose(4)
    fprintf('Reading binoc took %.2f, %d attempts\n',DATA.readdur,n);
end
if readflag
    DATA = GetState(DATA,'OpenPipe');
end
SetGui(DATA,'set');

    
function DATA = GetState(DATA, caller, verbose)
    if nargin < 3
        verbose = 0;
    end
    if nargin < 2
        caller = 'unknown';
    end
    if DATA.network
        str = [DATA.ip 'getstate'];
        ts = now;
        myprintf(DATA.frombinocfid,'Getstate caller: %s\n,',caller)
        [bstr, status] = urlread(str);
        a = mytoc(ts); %getting here is fast. Its interpretline that is slow.
        if isempty(bstr)
            vergwarning('GetState Return is Emtpy');
        elseif strncmp('STATEX',bstr,6)
            vergwarning('GetState Fault Rectified');
        end
         [DATA,~,~,nlines] = InterpretLine(DATA, bstr,'fromgetstate');
         if verbose
             fprintf('Read/Interpret took %.3f,%.3f (%d lines)\n',a,mytoc(ts),nlines);
         end
         myprintf(DATA.frombinocfid,'Getstate took %.3f,%.3f\n,',a,mytoc(ts))
    else
        outprintf(DATA,'QueryState\n');
        tic; DATA = ReadFromBinoc(DATA);toc
    end

function DATA = SetTrial(DATA, T)
    Trial = T;
    T.Start = now;
    nt = DATA.nt;
    if length(T.sv) && isfield(DATA,'exptype')
    DATA.Trials(nt).(DATA.exptype{1}) = T.sv(1);
    end
    if length(T.sv) > 1 && length(DATA.exptype) > 1
    DATA.Trials(nt).(DATA.exptype{2}) = T.sv(2);
    end
    if length(T.sv) > 2 && length(DATA.exptype) > 2
    DATA.Trials(nt).(DATA.exptype{3}) = T.sv(3);
    end
    f = fields(T);
    for j = 1:length(f)
        DATA.Trials(nt).(f{j}) = T.(f{j});
    end
    
    if isfield(DATA.binoc{1},'id')
        if ~isfield(T,'id') | isempty(T.id) | DATA.binoc{1}.id > T.id
        DATA.Trials(nt).id = DATA.binoc{1}.id;
        end
    end
    f = {'se'};
    for j = 1:length(f)
        if isfield(DATA.binoc{1},f{j})
            DATA.Trials(nt).(f{j}) = DATA.binoc{1}.(f{j});
        end
    end
    if ~isfield(DATA.Trials,'Nf') || isempty(DATA.Trials(nt).Nf)
        DATA.Trials(nt).Nf = NaN;
    end
    if ~isfield(DATA.Trials,'nf') || isempty(DATA.Trials(nt).nf)
        DATA.Trials(nt).nf = NaN;
    end
    if ~isfield(DATA.Trials,'dur') || isempty(DATA.Trials(nt).dur)
        DATA.Trials(nt).dur = NaN;
    end
    
    
function val = ReadVal(s, DATA)

    ccodes = {'ro' 'rx' 'ry'};
    truecodes = {'Ro' 'Rx' 'Ry'};
    
    id = strfind(s,'=');
    if isempty(id)
        s = s(3:end);
    else
        s = s(id(1)+1:end);
    end
    [val, n] = sscanf(s,'%f');
    if n == 0
        val = sscanf(s(3:end),'%f');
        if isempty(val)
            val = 0;
        end
        if sum(strcmp(s(1:2),ccodes))
            id = find(strcmp(s(1:2),ccodes));
            val = val + DATA.binoc{1}.(truecodes{id});
        elseif find(strcmp(s(1:2),{'Rx' 'Ry'}))
            val = val + DATA.binoc{1}.(s(1:2));
        end
    end

        
function S = SetField(S, F, x)
%Sets a field only if is doesn't exist
if ~isfield(S,F)
    S.(F) = x;
else
    x = S.(F);
end
    
function DATA = SetDefaults(DATA)

scrsz = get(0,'Screensize');
DATA = SetField(DATA,'ip','http://localhost:1110/');
DATA.lastreadtime = now;
DATA.smrfile = '';
DATA.showstatus.update = 1;
DATA.network = 2;
DATA.lastmsg = '';
DATA.errors(1) = 0; %keep track of  erros received, so only acknowlge first
DATA.confused = 0;
DATA.servodata.alldepths = [];
DATA.servodata.alltimes = [];
DATA.servodata.stepsize = 0;
DATA.Trial.rw = 0;
DATA.Comments = [];
DATA.state.stimfileerrs = 0;
DATA.state.stimfile = '';
DATA.state.query = 0;
DATA.state.dlgup = 0;
DATA.state.penwarned = 0;
DATA.currentrw = 0;
DATA.Header.dailyvars = {'id' 'se' 'ed' 'Rx' 'Ry' 'Ro' 'Rw' 'Rh' 'Xp' 'Yp' 'Pn' 'pe' 'Electrode' 'hemi'...
                        'ui' 'ePr' 'eZ' 'monkey' 'coarsemm' 'adapter' 'Trw' 'Tg' 'nT' 'Tb' 'uf' 'fx' 'fy' 'so' 'oldrf'};


DATA.newbinoc = 2;
DATA.ready = 0;
DATA.timerperiod = 0.05;
DATA.runsequence = 0;
DATA.canceltimer = 0;
DATA.pausereading = 0;  %stop timer driven reads when want to control
DATA.binocisup = 0;
DATA.pcdrive = '';
DATA.savestrs = 0;
DATA.restartbinoc = 0;
DATA.vergversion=vergversion();
DATA.matlabwasrun=0;
DATA.matexpres = [];
DATA.plotexpts = [];
DATA.completions = {};
DATA.helpstrs = {};
DATA.newchar = 0;
DATA.newexptdef = 1;  %when load 1st file, don't do any resets
DATA.netmatdir= [GetFilePath('binoclean') '/matlab'];
DATA.localmatdir='/local/matlab';
DATA.pausetime = 0;
DATA.readpause = 0;
DATA.rptexpts = 0;
DATA.verg.layoutfile = '/local/verg.layout';
DATA.font.FontSize = 14;
DATA.font.FontName = 'Arial';
DATA.Coil.gain= [0 0 0 0];
DATA.Coil.phase= [4 4 4 4];
DATA.Coil.offset= [0 0 0 0];
DATA.Coil.so = [0 0 0 0];
DATA.Coil.CriticalWeight = 0;
DATA.Coil.Weight = 0;
DATA.Coil.lastwt = 0;
DATA.layoutfile = '/local/verg.layout';
DATA.exptnextline = 0;
DATA.exptstoppedbyuser = 0;
DATA.pipelog = 0;
DATA.seqline = 0;
DATA.listmodified = [0 0 0];
DATA.Trial.Trial = 1;
DATA.windowcolor = [0.8 0.8 0.8];
DATA.Trial.sv = [];
DATA.psych.show = 1;
DATA.psych.trialresult = 0;
DATA.psych.showblocks = 0;
DATA.psych.blockmode = 'All';
DATA.psych.crosshairs = 1;
DATA.psych.collapse(2) = 0;
DATA.psych.collapse(3) = 1;
DATA.psych.blockid = [];
DATA.silentoption.NotSet = 1;
DATA.overcmds = {};
DATA.exptstimlist = { {} {} {} };
DATA.customstimlist = [0 0 0];
DATA.stimtypenames = {'fore' 'back' 'ChoiceU' 'ChoiceD'};
DATA.Statuslines = {};
DATA.statusitem = -1;
DATA.nt = 1;
DATA.exptype = [];
DATA.nexpts = 0;
DATA.Expts = {};
DATA.Trials = [];
DATA.showxy = [1 1 1]; %XY R, L, Conjugate crosses
DATA.currentstim = 1;  %foregr/backgre/Choice Targest
DATA.xyfsdvals = [1 2 5 10 20 40];
DATA.optionflags.ts = 0;
DATA.optionstrings.ts = 'Wurtz Task';
DATA.optionflags.do = 0;
DATA.optionstrings.do = 'Go';
DATA.showflags.ts = 1;
DATA.showflags.cf = 1;
DATA.showflagseq{1} = 'ts';
DATA.showflagseq{2} = 'cf';
DATA.stimflags{1}.pc = 1;
DATA.stimflags{1}.nc = 1;
DATA.stimflagnames.nc = 'Black Dots';
DATA.stimflagnames.pc = 'White Dots';
if ~isfield(DATA,'verbose')
    DATA.verbose = [0 0 0 0 0 1 0];
end
DATA = SetField(DATA,'matexpt','');
DATA = SetField(DATA,'perfmonitor',0);
DATA = SetField(DATA,'togglecodesreceived',0);
DATA = SetField(DATA,'autoreopen', 0);
DATA = SetField(DATA,'autorestart', 0);
DATA = SetField(DATA,'draintimeout', 4);
DATA = SetField(DATA,'pausetimeout', 30);
DATA = SetField(DATA,'nowarning',0);

DATA.commands = {};
DATA.commandctr = 1;
DATA.commandlines = [];
DATA.historyctr = 0;
DATA.inexpt = 0;
DATA.binoc{1}.uf = '';
DATA.electrodestrings = {'NotSet'};
DATA.userstrings = {'bgc' 'ali' 'ink' 'agb' 'pla' 'sid'};
DATA.monkeystrings = {'Icarus' 'Junior Barnes' 'Lemieux' 'Pepper' 'Rufus' };
DATA.monkeystrs = {'ica' 'jbe' 'lem' 'ppr' 'ruf' };
DATA.humanstrs = {'bgc' 'sid' 'pla' 'ink' 'agb'};

%make sure some fields exist. Be careful with this - if verg is resstarte
%and binoc is running, these might override what is in binoc
DATA.binoc{1}.Electrode = '';
DATA.binoc{1}.monkey = '';
DATA.binoc{1}.lo = '';
DATA.binoc{1}.nt = 1;
DATA.binoc{1}.st = ''; % make sure this field comes early
DATA.penid = 0;
DATA.stimulusnames{1} = 'none';
DATA.stimulusnames{4} = 'grating';
DATA.stimulusnames{3} = 'rds';
DATA.stimulusnames{2} = 'gabor';
DATA.stimulusnames{5} = 'bar';
DATA.stimulusnames{6} = 'circle';
DATA.stimulusnames{7} = 'rectangle';
DATA.stimulusnames{8} = 'test';
DATA.stimulusnames{9} = 'square';
DATA.stimulusnames{10} = 'probe';
DATA.stimulusnames{11} = '2grating';
DATA.stimulusnames{12} = 'cylinder';
DATA.stimulusnames{13} = 'corrug';
DATA.stimulusnames{14} = 'sqcorrug';
DATA.stimulusnames{15} = 'twobar';
DATA.stimulusnames{16} = 'rls';
DATA.redundantcodes = {'Bh' 'Bc' 'Bs' 'Op' 'Pp' 'sO' 'bO' 'aOp' 'aPp' 'O2' 'lf' 'rf'};
DATA.stimlabels = {'Foregnd' 'Backgnd' 'ChoiceU/R' 'ChoiceD/L'};

DATA.badnames = {'2a' '4a' '72'};
DATA.badreplacenames = {'afc' 'fc4' 'gone'};

DATA.comcodes = [];

%window names must be at least 8 chars
DATA.windownames = {'vergwindow' 'optionwindow' 'softoffwindow'  'codelistwindow' 'statuswindow' 'logwindow' 'helpwindow' 'sequencewindow' 'penlogwindow' 'electrodewindow' 'commentwindow' 'showpenwin'};
DATA.vergonlycodes = {'timerperiod' 'autoreopen' 'pausetimeout'};
DATA.winpos{1} = [10 scrsz(4)-480 300 450];
DATA.winpos{2} = [10 scrsz(4)-680 400 50];  %options popup
DATA.winpos{3} = [600 scrsz(4)-100 600 150]; %softoff
DATA.winpos{4} = [600 scrsz(4)-600 400 500]; %code list
DATA.winpos{5} = [600 scrsz(4)-100 500 400]; %status
DATA.winpos{6} = [600 scrsz(4)-100 400 100]; %log
DATA.winpos{7} = [600 scrsz(4)-100 400 100]; %help
DATA.winpos{8} = [600 scrsz(4)-100 400 100]; %sequence
DATA.winpos{9} = [600 scrsz(4)-100 400 100]; %Penetraation log
DATA.winpos{10} = [600 scrsz(4)-100 400 100]; %Electrode Moving
DATA.winpos{11} = [600 scrsz(4)-100 400 100]; %Electrode Moving
DATA.winpos{12} = [600 scrsz(4)-100 400 100]; %Electrode Moving

DATA.outid = 0;
DATA.inid = 0;
DATA.cmdfid = 0;
DATA = SetField(DATA,'frombinocfid', 0);
DATA.tobinocfid = 0;
DATA.incr = [0 0 0];
DATA.nstim = [1 1 1];
DATA.quickexpts = [];
DATA.helpfiles = [];
DATA.stepsize = [20 10];
DATA.stepperpos = -2000;
DATA.tag.stepper = 'Stepper';
DATA.tag.softoff = 'Softoff';
DATA.tag.options = 'Options';
DATA.tag.penlog = 'Penetration Log';
DATA.tag.monkeylog = 'Monkey Log';
DATA.tag.codes = 'Codelist';
DATA.tag.psych = 'VergPsych';
DATA.tag.status = 'StatusWindow';
DATA.tag.plotwin = 'PlotWindow';
%make sure the size field is defined, then it comes before wi,hi
DATA.comcodes(1).label = 'Size';
DATA.comcodes(1).code = 'sz';
DATA.comcodes(1).const = 1;
DATA.comcodes(1).type = 'num';
DATA.comcodes(1).group = 1;
DATA.strcodes(1).label = 'Old Monitor file';
DATA.strcodes(1).code = 'oldmonitor';
DATA.strcodes(1).icode = 0; 

%Need to initialize values that verg might try to use 
%before reading values from binoclean
%Also define fields here that need to be send early in stream
% like sz - send before wi/hi so that is doen not override
DATA.binoc{1}.xo = 0;
DATA.binoc{1}.sz = 0;
DATA.binoc{2}.xo = 0;
DATA.binoc{1}.ePr = 0;
DATA.binoc{1}.adapter = 'None';
DATA.binoc{1}.uf = '';
DATA.binoc{1}.ei = '0';
DATA.binoc{1}.i2 = '0';
DATA.binoc{1}.i3 = '0';
DATA.binoc{1}.nr = 1;
DATA.binoc{1}.rw = 0;
DATA.binoc{1}.seqpause = 10;
DATA.binoc{1}.ereset = 'NotSet';
DATA.binoc{1} = SetField(DATA.binoc{1},'magic',0);

DATA.binocstr.monitor = '/local/monitors/Default';
DATA.completestr = '';
DATA.expts{1} = [];
DATA.expts{2} = [];
DATA.expts{3} = [];
DATA.expmenuvals{1} = [];
DATA.expmenuvals{2} = [];
DATA.expmenuvals{3} = [];
DATA.extypes{1} = [];
DATA.extypes{2} = [];
DATA.extypes{3} = [];
DATA.exptype{1} = 'e0';
DATA.exptype{2} = 'e0';
DATA.exptype{3} = 'e0';
DATA.stimtype(1) = 1;
DATA.stimtype(2) = 1;
DATA.nextras = [0 0 0];

DATA.electrodeid = 1;
DATA.mean = [0 0 0];
DATA.incr = [0 0 0];
DATA.servoport = []; %name of tty2 in /local/binoc.setup
DATA.servofig = 0;  %Figre # of servo control window.
DATA = ReadSetupFile(DATA, '/local/binoc.setup');
DATA = ReadStimFile(DATA, '/local/verg.setup');
[DATA.helpstrs, DATA.helpkeys] = ReadHelp(DATA);

if ~isfield(DATA,'matexpvars') %may have been set in verg.setup
    DATA.matexpvars = [];
end

for j = 1:3 
    DATA.expmenucodes{j} = {};
    DATA.expstrs{j} = {};
    DATA.expmenuvals{j} = [];
end
for j = 1:length(DATA.comcodes)
    if ismember(DATA.comcodes(j).const,DATA.extypes{1})
        DATA.expstrs{1} = {DATA.expstrs{1}{:} DATA.comcodes(j).label};
        DATA.expmenuvals{1} = [DATA.expmenuvals{1} DATA.comcodes(j).const];
        DATA.expmenucodes{1} = {DATA.expmenucodes{1}{:} DATA.comcodes(j).code};
    end
    if ismember(DATA.comcodes(j).const,DATA.extypes{2})
        DATA.expstrs{2} = {DATA.expstrs{2}{:} DATA.comcodes(j).label};
        DATA.expmenuvals{2} = [DATA.expmenuvals{2} DATA.comcodes(j).const];
    end
    if ismember(DATA.comcodes(j).const,DATA.extypes{3})
        DATA.expstrs{3} = {DATA.expstrs{3}{:} DATA.comcodes(j).label};
        DATA.expmenuvals{3} = [DATA.expmenuvals{3} DATA.comcodes(j).const];
    end
end
txt = scanlines('/local/verg.setup');
id = find(strncmp('autostart',txt,9));
if ~isempty(id) && DATA.autorestart == 0
    DATA.autorestart = 2;
end

function [strs, Keys] = ReadHelp(DATA)
 
    strs = {};
    Keys = [];
    helpfile = [DATA.localmatdir '/helpstrings.txt'];
    fid = fopen(helpfile,'r');
    Keys.extras =[];
    Keys.cmdcode = [];
    if fid > 0
    a = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    txt = a{1};
    code = [];
    lastcode = code;
%In helpstrings.txt, lines beginning with
%# provide additional help on the preceding code. Displayed when line is selectd
%+ help on option codes
%' ignored
    for j = 1:length(txt)
        code = regexprep(txt{j},'\s.*','');
        if isempty(code) || code(1) == ''''
        elseif code(1) == '+'  %optionflag help
            str = regexprep(txt{j},code,'');
            code = code(2:end);
            Keys.options.(code) = str;
            lastcode = code;
        elseif txt{j}(1) == '#' 
            if ~isempty(lastcode)
                if isfield(Keys.extras,lastcode)
                    Keys.extras.(lastcode){end+1} = txt{j}(2:end);
                else
                    Keys.extras.(lastcode){1} = txt{j}(2:end);
                end
            else
                code
            end
        elseif ~isempty(code) && txt{j}(1) ~= '#'
            if isfield(strs,code) && DATA.verbose(4)
                fprintf('%s help duplicated\n',code);
            end
            str = regexprep(txt{j},code,'');
            if code(1) == '!'
                code = code(2:end);
                Keys.cmdcode.(code) = '!';
            end
            if strfind(str,'#')
                keyword = regexprep(str,'.*#','');
                Keys.KeyWords.(code) = keyword;
            else
                Keys.KeyWords.(code) = '';
            end     
            if strfind(str,'\#')
                str = strrep(str,'\#','#');
            else
                str = regexprep(str,'#.*','');
            end
            if isfield(Keys.cmdcode,code)
                strs.(code) =  [Keys.cmdcode.(code) str];
            else
                strs.(code) =  str;
            end
            lastcode = code;
        end
    end
    end
    
    
 function DATA = SetExptMenus(DATA)

    
    
if isfield(DATA,'exptstimlist')
    for j = 1:length(DATA.exptstimlist)
        S = DATA.exptstimlist{j};
        for m = 1:length(S)
            [val,n] = sscanf(S{m},'%f');
            if n == 1
            DATA.expvals{j}(m) = val;
            end
        end
    end
end



for m = 1:3
    if length(DATA.expmenuvals{m}) > 0
        DATA.expstrs{m} = {DATA.comcodes(DATA.expmenuvals{m}).label};
        DATA.expmenucodes{m} = {DATA.comcodes(DATA.expmenuvals{m}).code};
    else
    DATA.expmenucodes{m} = {};
    DATA.expstrs{m} = {};
    DATA.expmenuvals{m} = [];
    for j = 1:length(DATA.extypes{m})
        id = find([DATA.comcodes.const] == DATA.extypes{m}(j));
        if length(id) == 1
            DATA.expstrs{m} = {DATA.expstrs{m}{:} DATA.comcodes(id).label};
            DATA.expmenuvals{m} = [DATA.expmenuvals{m} DATA.comcodes(id).const];
            DATA.expmenucodes{m} = {DATA.expmenucodes{m}{:} DATA.comcodes(id).code};
        end
    end
    end
end

   
if isfield(DATA,'toplevel') %GUI is up
it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt1List');
set(it,'string',DATA.expstrs{1});
it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt2List');
set(it,'string',DATA.expstrs{2});
it = findobj(allchild(DATA.toplevel),'flat','Tag','Expt3List');
e2 = get(it,'value');
if e2 > length(DATA.expstrs{3})
    set(it,'value',1);
end
set(it,'string',DATA.expstrs{3});
end

function ShowStatus(DATA)

    persistent oldcount;
    if isempty(oldcount)
        oldcount = 0;
    end
    status = 1;
    if DATA.inexpt  && DATA.nexpts > 0
        str = datestr(DATA.Expts{DATA.nexpts}.Start);
        str = ['Started ' str(13:17)];
    elseif DATA.nexpts > 0
        if isfield(DATA.Expts{DATA.nexpts},'End')
            str = datestr(DATA.Expts{DATA.nexpts}.End);
            str = ['Ended ' str(13:17)];
        else
            str = sprintf('Binoc seems not to be in Expt, but verg didnt get end of Expt %d',DATA.nexpts);
        end
    else
        status = 0;
        str = ['No Trials Run'];
    end
    if isfield(DATA,'trialcounts') && length(DATA.trialcounts) > 7
    s = sprintf('Trials %d/%d Bad%d Late%d  Rw%.1f (%.2f) Ex:%d/%d %s',...
    DATA.trialcounts(1),DATA.trialcounts(2),DATA.trialcounts(3),DATA.trialcounts(4),...
    DATA.trialcounts(8),DATA.Trial.rw,DATA.trialcounts(6),DATA.trialcounts(7),str);
    else
        s = str;
    end

if isfield(DATA,'toplevel')
    set(DATA.toplevel,'Name',s);
end
%trialcounts(2) is total trials according to binoc
%print out trails counts when this increases.
if isfield(DATA,'trialcounts')
    if (status > 0 && DATA.verbose(6) && DATA.trialcounts(2) > oldcount) || DATA.verbose(3)
        oldcount = DATA.trialcounts(2);
        fprintf('%s\n',s);
    end
else
    if DATA.verbose(4) %verg state
        cprintf('red','No Trialcounts\n');
    end
end


function DATA = InitInterface(DATA)

    scrsz = get(0,'Screensize');
    cntrl_box = figure('Position', DATA.winpos{1},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{1},'Name',DATA.name,'menubar','none');
                set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
                set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);

%?? do we still need to do this? 
    if isfield(DATA.showflags,'do') && DATA.togglecodesreceived == 0 
    DATA.showflags = rmfield(DATA.showflags,'do');
    end
    f = fields(DATA.showflags);
    nc = 5;
    tagc = 6;%# of columns for toggles
    nr = 19 + ceil((1+length(f))./tagc);
    cw = 0.99/nc;
    DATA.toplevel = cntrl_box;
    
    lst = jcontrol(gcf,'javax.swing.JTextField',...
                    'Units','normalized',...
                    'Position',[0.01 0.01 0.98 1./nr]);

    set(lst, 'KeyPressedCallback', {@jTextKey});
    lst.setFocusable(true);
    lst.setFocusTraversalKeysEnabled(false)
%    lst.KeyPressedCallback = @jTextKey
    lst.setBackground(java.awt.Color(1,1,1))
%     set(lst, 'String', '',...            
%         'HorizontalAlignment','left',...
%         'Callback', {@TextEntered}, 'Tag','NextButton',...
%         'units','norm', 'Position',[0.01 0.01 0.98 1./nr],...
%          'KeyReleasedCallBack', @JtextKey);
    DATA.txtui = lst;
%        'KeyPressFcn',@TextKey,
    
    
    bp = [0.01 1.01./nr 3.5./nc 6/nr];
    lst = uicontrol(gcf, 'Style','list','String', 'Command History',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
        'Tag','CommandHistory',...
        'callback',@TextList,...
'units','norm', 'Position',bp);
   DATA.txtrec = lst;
   
    xp = bp(1)+bp(3);
    bp(1) = bp(1)+bp(3);
    bp(2) = 2./nr;
    bp(3) = cw/3;
    bp(4) = 1./nr;
    if isfield(DATA.binoc{1},'xyfsd')
        [a,j] = min(abs(DATA.binoc{1}.xyfsd - DATA.xyfsdvals));
    else
        j = 1;
    end
    uicontrol(gcf,'style','text','string','FSD',  'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3)=0.99-bp(1);
    uicontrol(gcf,'style','pop','string',num2str(DATA.xyfsdvals'), ...
        'units', 'norm', 'position',bp,'value',j,'Tag','FSD','callback',{@SetExpt, 'fsd'});

     bp(3) = cw/2;
    bp(2) = 3./nr;
    bp(1) = xp;
    ai = uicontrol(gcf,'style','checkbox','string','XYR',  'Tag', 'XYR', 'value', DATA.showxy(1), 'units', 'norm', 'position',bp, 'callback', {@OtherToggles, 'XYR'});
    
    bp(1) = bp(1)+bp(3);
    uicontrol(gcf,'style','checkbox','string','XYL',  'Tag', 'XYL', 'value', DATA.showxy(2), 'units', 'norm', 'position',bp, 'callback', {@OtherToggles, 'XYL'});
    bp(1) = bp(1)+bp(3);
    uicontrol(gcf,'style','checkbox','string','XYB',  'Tag', 'XYB', 'value', DATA.showxy(3), 'units', 'norm', 'position',bp, 'callback', {@OtherToggles, 'XYB'});

    if GetOption(DATA,'py')
        bp(1) = bp(1)+bp(3);
        bp(3) = 0.25;
        id = find(strcmp(DATA.binoc{1}.monkey,DATA.monkeystrs));
        if isempty(id)
            id = 1;
        end
        h = uicontrol(gcf,'style','pop','string',DATA.humanstrs, ...
            'units', 'norm', 'position',bp,'value',id,'Tag','Monkey','callback',{@MenuGui});                
        gui.place(h,'up', ai);
        bp(3) = 0.9./nc;
        a = uicontrol(gcf,'style','text','string','Subject', ...
            'units', 'norm', 'position',bp);
        gui.place(a, 'up',h);
    end
    
    bp(1) = 0.01;
    bp(2) = 8./nr;
    bp(3) = 0.1;
    bp(4) = 1./nr;
    uicontrol(gcf,'style','text','string','File',  'units', 'norm', 'position',bp);

    
    
    bp(1) = bp(1)+bp(3);
    bp(3) = 1-bp(1)-0.1;
    uicontrol(gcf,'style','edit','string',DATA.binoc{1}.uf, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','DataFileName','callback',{@TextGui, 'uf'});
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.1;
            uicontrol(gcf,'style','pushbutton','string','Open', ...
        'Callback', {@OpenUffFile, 1}, 'Tag','UffButton',...
        'units', 'norm', 'position',bp,'value',1);
    
    
    bp(1) = 0.01;
    bp(2) = 7./nr;
    bp(3) = 0.15;
    bp(4) = 1./nr;
    uicontrol(gcf,'style','text','string','Comment',  'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 1-bp(1);
    uicontrol(gcf,'style','edit','string','', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Comment','callback',{@TextGui, 'cm'});

    
    bp(2) = 11./nr;
    bp(1) = 0.01;
    bp(4) = 1./nr;
    bp(3) = cw;
        uicontrol(gcf,'style','pushbutton','string','Finish', ...
        'Callback', {@RunButton, 2}, 'Tag','StopButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(2) = 10./nr;
    bp(3) = cw/2;
         uicontrol(gcf,'style','Text','string','Rpt', ...
        'units', 'norm', 'position',bp,'value',1);
    bp(2) = 10./nr;
    bp(1) = cw/2;
         uicontrol(gcf,'style','edit','string',num2str(DATA.rptexpts), ...
             'Tag','RptExpts','Callback',@TextGui,...
        'units', 'norm', 'position',bp,'value',1);

    bp(2) = 12./nr;
    bp(1) = 0.01;
    bp(4) = 1./nr;
    bp(3) = cw/4;
         uicontrol(gcf,'style','Text','string','N', ...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(3);
    bp(3) = 0.75*cw;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.nr), 'Tag', 'binoc.nr', 'units', 'norm',...
        'callback',{@TextGui, 'nr'},'position',bp);
    
    
    
    cmenu = uicontextmenu;
    uimenu(cmenu,'label','Apply to Expt 1 (cntrl-g)','Callback',{@EditValsMenu, 'apply'});
    uimenu(cmenu,'label','Cancel (Esc)','Callback',{@EditValsMenu, 'cancel'});
    
    bp(3) = 1./nc;
    bp(2) = 9./nr;
    bp(4) = 4/nr;
    bp(1) = 0.02+cw;
 %elsewhere use findobj to get this item, so only use the tag once
    a= uicontrol(gcf,'style','edit','string',DATA.exptstimlist{1}, 'min',1,'max',5,...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1StimList','keypressfcn',@EditText,'callback',@TextCallback);
    set(cmenu,'UserData',a,'tag','Expt1StimList');
    set(a,'uicontextmenu',cmenu);

    
    cmenu = uicontextmenu;
    uimenu(cmenu,'label','Apply to Expt2 (cntrl-g)','Callback',{@EditValsMenu, 'apply2'});
    uimenu(cmenu,'label','Cancel (Esc)','Callback',{@EditValsMenu, 'cancel'});
    bp(1) = bp(1)+bp(3)+0.01;
    a = uicontrol(gcf,'style','edit','string',num2str(DATA.nstim(2)),  'min',1,'max',5,...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2StimList','keypressfcn',@EditText);
    bp(1) = bp(1)+bp(3)+0.01;
    set(cmenu,'UserData',a,'tag','Expt2StimList');
    set(a,'uicontextmenu',cmenu);

    cmenu = uicontextmenu;
    uimenu(cmenu,'label','Apply to Expt3 (cntrl-g)','Callback',{@EditValsMenu, 'apply3'});
    uimenu(cmenu,'label','Cancel (Esc)','Callback',{@EditValsMenu, 'cancel'});
    a = uicontrol(gcf,'style','edit','string',num2str(DATA.nstim(3)),  'min',1,'max',5, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3StimList','keypressfcn',@EditText,'callback',@TextCallback);
    set(cmenu,'UserData',a,'tag','Expt3StimList');
    set(a,'uicontextmenu',cmenu);
    
    
    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = cw;
    bp(4) = 1.2/nr;
    uicontrol(gcf,'style','text','string','N stim',  'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = cw;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.nt), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1Nstim','callback',{@TextGui, 'nt'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.nstim(2)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2Nstim','callback',{@TextGui, 'n2'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.nstim(3)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3Nstim','callback',{@TextGui, 'n3'});

    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = cw;
    uicontrol(gcf,'style','text','string','Incr',  'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = cw;
    uicontrol(gcf,'style','edit','string',DATA.binoc{1}.ei, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1Incr','callback',{@TextGui, 'ei'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',DATA.binoc{1}.i2, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2Incr','callback',{@TextGui, 'i2'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',DATA.binoc{1}.i3, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3Incr','callback',{@TextGui, 'i3'});

    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = cw;
    uicontrol(gcf,'style','text','string','Mean',  'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = cw;
    uicontrol(gcf,'style','edit','string',num2str(DATA.mean(1)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1Mean','callback',{@TextGui, 'em'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.mean(2)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2Mean','callback',{@TextGui, 'm2'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.mean(3)), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3Mean','callback',{@TextGui, 'm3'});

    
    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = cw;
    
    uicontrol(gcf,'style','pushbutton','string','Run', ...
        'Callback', {@RunButton, 1}, 'Tag','RunButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = cw;
    uicontrol(gcf,'style','pop','string',DATA.expstrs{1}, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt1List','callback',{@SetExpt, 'et'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pop','string',DATA.expstrs{2}, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt2List','callback',{@SetExpt, 'e2'});
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pop','string',DATA.expstrs{3}, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Expt3List','callback',{@SetExpt, 'e3'});
    
    bp(1) = 0.01;
    bp(2) = bp(2)+bp(4);
    bp(3) = cw;
    a = uicontrol(gcf,'style','text','string','Foreground',  'units', 'norm', 'Tag','CurrentStimLabel','position',bp);     
    cmenu = uicontextmenu;
    uimenu(cmenu,'label','Edit ForeGround','Callback',{@StimulusSelectMenu, 'fore'});
    uimenu(cmenu,'label','Edit BackGround','Callback',{@StimulusSelectMenu, 'back'});
    uimenu(cmenu,'label','Edit Choice Icon (U/R)','Callback',{@StimulusSelectMenu, 'ChoiceU'});
    uimenu(cmenu,'label','Edit Choice Icon (D/L)','Callback',{@StimulusSelectMenu, 'ChoiceD'});
    set(a,'uicontextmenu',cmenu);
    
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = cw;
    uicontrol(gcf,'style','pop','string',DATA.stimulusnames, ...
        'units', 'norm', 'position',bp,'value',DATA.stimtype(1),'Tag','ForegroundType','callback',{@SetExpt, 'st'});
    bp(1) = bp(1)+bp(3)+0.01;
    a = uicontrol(gcf,'style','text','string','Backgnd',  'units', 'norm', 'Tag','BackgrStimLabel','position',bp);     
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pop','string',DATA.stimulusnames, ...
        'units', 'norm', 'position',bp,'value',DATA.stimtype(2),'Tag','BackgroundType','callback',{@SetExpt, 'bs'});

    

    
        
    cmenu = uicontextmenu;
    f = fields(DATA.optionstrings);
    of = fields(DATA.optionflags);
    amenu = uimenu(cmenu,'Label','General');
    for j = 1:length(f)
        if strncmp(of{j},'lbl',3)
            amenu = uimenu(cmenu,'Label',DATA.optionstrings.(f{j}));
        else
            h = uimenu(amenu,'label',DATA.optionstrings.(f{j}),'Callback',{@OptionMenu, f{j}},'tag',f{j});
            if DATA.optionflags.(f{j})
                set(h,'checked','on');
            end
        end
    end
    set(cmenu,'tag','OptionContextMenu');

    
    
    bp(1) = 0.01;
    bp(2) = 1-1/nr;
    bp(4) = 1./nr;
    bp(3) = 1./tagc;
    uicontrol(gcf,'style','checkbox','string','Go', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','do','callback',@GoToggle);
    f = fields(DATA.showflags);
    allf = fields(DATA.optionflags);
    if isfield(DATA,'showflagseq')
        f = DATA.showflagseq;
        nrows = ceil(length(f)./tagc)./nr; %fraction of window height needed for tags
    else
        nrows = 2/nr;
    end
    
    
    f = f(find(~strcmp('do',f)));
    ymin = 1-1./nr - nrows;
    ymin = 1 - nrows;
    for j = 1:length(f)
        id = find(strcmp(f{j},allf));
        if length(id) == 1
            str = DATA.optionstrings.(allf{id});
        else
            str = num2str(j);
        end
        bp(2) = bp(2)-bp(4);
        if bp(2) < ymin
            bp(1) = bp(1)+bp(3);
            bp(2) = 1- 1./nr;
        end
        if isfield(DATA.optionflags,f{j})
            h = uicontrol(gcf,'style','checkbox','string',str, ...
                'units', 'norm', 'position',bp,'value',DATA.optionflags.(f{j}),'Tag',f{j},'callback',{@HitToggle, f{j}});
            set(h,'uicontextmenu',cmenu);
        end

    end
    bp(3) = 1/nc;
    
    hm = uimenu(cntrl_box,'Label','File','Tag','BinocFileMenu');
    uimenu(hm,'Label','Close','Callback',{@verg, 'close'});
    uimenu(hm,'Label','Close Verg and Binoc','Callback',{@MenuHit, 'bothclose'});
    uimenu(hm,'Label','Restart Binoc','Callback',{@MenuHit, 'restartbinoc'});
    uimenu(hm,'Label','Copy Logs to PC/Network','Callback',{@MenuHit, 'copylogs'});
    uimenu(hm,'Label','Save','Callback',{@SaveFile, 'current'});
    uimenu(hm,'Label','Save As...','Callback',{@SaveFile, 'saveas'});
    if DATA.state.stimfileerrs
        sname = strrep(DATA.state.stimfile,fileparts(DATA.state.stimfile),'');
        sname = sname(2:end);
        uimenu(hm,'Label',sprintf('Fix %s',sname),'Callback',{@SaveFile, 'fix'});
    end
    sm = uimenu(hm,'Label','Preferences');
    uimenu(sm,'Label','Save Layout','Callback',{@SaveFile, 'layout'});
    uimenu(sm,'Label','Choose Font','Callback',{@MenuHit, 'choosefont'});
    uimenu(sm,'Label','Update local .m files','Callback',{@SaveFile, 'update'});
    uimenu(sm,'Label','Push local .m files to network binoclean','Callback',{@SaveFile, 'push'});
    sm = uimenu(hm,'Label','Today Slots');
    for j = 1:8
    uimenu(sm,'Label',['Expt' num2str(j)],'Callback',{@SaveSlot, j});
    end
    sm = uimenu(hm,'Label','Recover','Callback',{@RecoverFile, 'toplist'});
    x = uimenu(sm,'Label','List','Callback',{@RecoverFile, 'list'});
    uimenu(sm,'Label','eo.stm','Callback',{@RecoverFile, 'eo.stm'});
    uimenu(sm,'Label','eb.stm','Callback',{@RecoverFile, 'eb.stm'});
    uimenu(sm,'Label','0.stm','Callback',{@RecoverFile, '0stim'});
    uimenu(sm,'Label','1.stm','Callback',{@RecoverFile, '2stim'});
    uimenu(sm,'Label','Just Read id/se from last','Callback',{@RecoverFile, 'loadlast'});
    set(DATA.toplevel,'UserData',DATA);
            RecoverFile(x,[],'list');    
    hm = uimenu(cntrl_box,'Label','Quick','Tag','QuickMenu');
    BuildQuickMenu(DATA, hm);
%    uimenu(hm,'Label','Stepper','Callback',{@StepperPopup});
%currently can do everything in binoc. Stick with this til need something
%new....

    
    hm = uimenu(cntrl_box,'Label','&Pop','Tag','PopMenu');
    subm = uimenu(hm,'Label','&Software Offset');
    uimenu(subm,'Label','&Null','Callback',{@SendStr, 'sonull'},'accelerator','E');
    uimenu(subm,'Label','Edit','Callback',{@SoftoffPopup, 'popup'});
    uimenu(subm,'Label','Clear','Callback',{@SendStr, '\clearsoftoff'});
    uimenu(hm,'Label','Run Sequence of expts','Callback',{@SequencePopup, 'popup'});
    uimenu(hm,'Label','Pause Expt','Callback',{@SendStr, '\pauseexpt'});
    uimenu(hm,'Label','Center stimulus','Callback',{@SendStr, 'centerstim'});
    uimenu(hm,'Label','Clear Drawn Lines','Callback',{@SendStr, '!clearlines'});
    uimenu(hm,'Label','Remember This RF','Callback',{@SendStr, '!saverf'});
    AddMarkMenu(uimenu(hm,'Label','Mark'));
    
    uimenu(hm,'Label','BlackScreen (shake)','Callback',{@MenuHit, 'setshake'},'accelerator','B');
    uimenu(hm,'Label','pipelog','Callback',{@MenuHit, 'pipelog'});
    uimenu(hm,'Label','Update Network Psych Files','Callback',{@MenuHit, 'updatepsych'});
    uimenu(hm,'Label','freereward','Callback',{@MenuHit, 'freereward'},'accelerator','R');
    uimenu(hm,'Label','Set Spike Display (Spike2)','Callback',{@SendStr,'spike2:setspkv'});
    uimenu(hm,'Label','Run One Trial','Callback',{@MenuHit, 'onetrial'},'accelerator','1');

    
    
    
%    hm = uimenu(cntrl_box,'Label','Mark');

    subm = uimenu(cntrl_box,'Label','&Windows');
    uimenu(subm,'Label','&Options','Callback',{@OptionPopup},'accelerator','O');
    uimenu(subm,'Label','Electrode Control','Callback',{@ElectrodePopup, 'popup'});
    uimenu(subm,'Label','Comments','Callback',{@CommentPopup, 'popup'});
    uimenu(subm,'Label','List Codes','Callback',{@CodesPopup, 'popup'},'accelerator','L');
    uimenu(subm,'Label','Monkey Log','Callback',{@MonkeyLogPopup, 'popup'});
    uimenu(subm,'Label','Penetration Log','Callback',{@PenLogPopup});
    uimenu(subm,'Label','Pen Log Contents','Callback',{@ShowPenLog, 'popup'});
    uimenu(subm,'Label','Psych Window','Callback',{@MenuHit, 'showpsych'});
    uimenu(subm,'Label','Status Lines','Callback',{@StatusPopup, 'popup'});
    uimenu(subm,'Label','Software Offset','Callback',{@SoftoffPopup, 'popup'});
    sm = uimenu(subm,'Label','Plots');
    uimenu(sm,'Label','Trial Durations Histogram (Day)','Callback',{@MenuHit, 'checkdurhist'});
    uimenu(sm,'Label','Trial Durations Histogram (Expt)','Callback',{@MenuHit, 'checkdurexpt'});
    uimenu(sm,'Label','Trial Duration sequence','Callback',{@MenuHit, 'checkdurseq'});


    subm = uimenu(cntrl_box,'Label','Pipes');
    uimenu(subm,'Label','Reopen Pipes','Callback',{@ReadIO, 6});
    uimenu(subm,'Label','reopenserial','Callback',{@SendStr, '\reopenserial'});
    uimenu(subm,'Label','Read','Callback',{@ReadIO, 1});
    uimenu(subm,'Label','GetState','Callback',{@ReadIO, 2});
    uimenu(subm,'Label','NewStart','Callback',{@ReadIO, 3});
    uimenu(subm,'Label','Stop Timer','Callback',{@ReadIO, 4});
    sm = uimenu(subm,'Label','Check Timer','Callback',{@CheckTimerHit, 0});
    sm = uimenu(subm,'Label','Start Timer','Callback',{@ReadIO, 5},'foregroundcolor',[0 0 0.5]);
    sm = uimenu(subm,'Label','Verbose');
    xm = uimenu(sm,'Label','To Binoc', 'Callback', {@SetVerbose, 1});
    SetMenuCheck(xm, DATA.verbose(1));
    xm = uimenu(sm,'Label','From Binoc','Callback', {@SetVerbose, 7});
    SetMenuCheck(xm, DATA.verbose(7));
    xm = uimenu(sm,'Label','Trial Data', 'Callback',{@SetVerbose, 3});
    SetMenuCheck(xm, DATA.verbose(3));
    xm = uimenu(sm,'Label','verg state', 'Callback',{@SetVerbose, 4});
    SetMenuCheck(xm, DATA.verbose(4));
    xm = uimenu(sm,'Label','IOTiming', 'Callback',{@SetVerbose, 5});
    SetMenuCheck(xm, DATA.verbose(5));
    xm = uimenu(sm,'Label','Trial Results', 'Callback',{@SetVerbose, 6});
    SetMenuCheck(xm, DATA.verbose(6));
    xm = uimenu(sm,'Label', 'Suppress All Popups', 'Callback', {@SetVerbose, 'nopopup'});
    SetMenuCheck(xm, DATA.nowarning);
    xm = uimenu(sm,'Label', 'Open Disk Log ToBinoc', 'Callback', {@SetVerbose, 'openlogto'});
    if DATA.tobinocfid > 0
        set(xm,'checked','on');
    end
    xm = uimenu(sm,'Label', 'Open Disk Log FromBinoc', 'Callback', {@SetVerbose, 'openlogfrom'});
    if DATA.frombinocfid > 0
        set(xm,'checked','on');
    end
    xm = uimenu(sm,'Label', 'All Off', 'Callback', {@SetVerbose, 0});
  
    sm = uimenu(subm,'Label','Try Pipes','Callback',{@ReadIO, 8},'foregroundcolor','r');
    sm = uimenu(subm,'Label','Restart Binoc between Expts','Callback',{@SetMenuToggle, 'restartbinoc'});
    sm = uimenu(subm,'Label','Send Expt Trigger to Spike2','Callback',{@SendStr, '!expttrigger'});
    sm = uimenu(subm,'Label','Send Spike Clear to Spike2','Callback',{@SendStr, 'cl'});
    uimenu(subm,'Label','Log Inputs','Callback',{@ReadIO, 'openlog'});
    sm = uimenu(subm,'Label','Tests');
    uimenu(sm,'Label','Run/Cancel','Callback',{@TestIO, 'cancel'});
    uimenu(sm,'Label','Freeze','Callback',{@TestIO, 'freeze'});
    uimenu(sm,'Label','Step','Callback',{@TestIO, 'step'});
    uimenu(sm,'Label','Save Image','Callback',{@TestIO, 'saveimage'});

    
    hm = uimenu(cntrl_box,'Label','Help','Tag','HelpMenu');
    DATA = AddHelpFiles(DATA); 
    BuildHelpMenu(DATA, hm);

    if DATA.network
        period = 0.005;
        period = DATA.timerperiod;  %don't need 10ms reads yet
    else
        period = 1;
    end
    DATA.timerobj = timer('timerfcn',{@CheckInput, DATA.toplevel},'Period',period,'executionmode','fixedspacing');
    setappdata(DATA.toplevel,'PauseReading',0);
    
    set(DATA.toplevel,'UserData',DATA);
    start(DATA.timerobj);

function ShowHelp(a,b,file)

  DATA = GetDataFromFig(a);
  cntrl_box = findobj('Tag',DATA.windownames{7},'type','figure');
  if isempty(cntrl_box)
  cntrl_box = figure('Position', DATA.winpos{7},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{7},'Name','HelpText','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
        set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
        set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);

    
    lst = uicontrol(gcf, 'Style','edit','String', 'HelpText',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
        'Tag','HelpText',...
'units','norm', 'Position',[0.01 0.01 0.99 0.99]);
  else
      lst = findobj(cntrl_box,'Tag','HelpText')
      figure(cntrl_box);
  end
  try
      fid = fopen(file,'r');
      fclose(fid)
      txt = textread(file,'%s','delimiter','\n');
      set(lst,'string',txt);
  catch
      fprintf('%s\n',lasterr)
  end
  
  function DATA = AddHelpFile(DATA, s)
    s = s(10:end);
        id = strfind(s,'"');
        if ~isempty(id)
            n = [];
            lbl = s(2:id(end)-1);
            if isfield(DATA.helpfiles,'label')
            n = find(strcmp(lbl,{DATA.helpfiles.label}));
            end
            if isempty(n)
                n = length(DATA.helpfiles)+1;
            end
            DATA.helpfiles(n).filename = s(id(end)+1:end);
            DATA.helpfiles(n).label = s(2:id(end)-1);
        end

  
  function DATA = AddHelpFiles(DATA, varargin)
      preifx = [];
      if isfield(DATA.binoc{1},'helpdir') && isdir(DATA.binoc{1}.helpdir)
          prefix = DATA.binoc{1}.helpdir;
      elseif isdir('/b/binoclean/help')
          prefix = '/b/binoclean/help';
      end
      if ~isempty(prefix)
          d = dir([prefix '/*.hlp']);
          for j = 1:length(d)
              if isempty(DATA.helpfiles) || ...
                      sum(cellstrfind({DATA.helpfiles.filename},d(j).name)) == 0 %new
                  if ~strncmp(d(j).name,'.#',2)
                      filename = [prefix '/' d(j).name];
                      DATA.helpfiles(end+1).filename = filename;
                      s = scanlines(filename);
                      if ~isempty(s) && length(s{1}) < 50
                          DATA.helpfiles(end).label = s{1};
                      else
                          DATA.helpfiles(end).label = d(j).name;
                      end
                  end
              end
          end
          
      end
  
   function BuildHelpMenu(DATA, hm)
        
      
   uimenu(hm,'Label','List All Codes','Callback',{@CodesPopup, 'popup'},'accelerator','L');
   uimenu(hm,'Label','List Codes with Help','Callback',{@CodesPopup, 'popuphelp'});
   for j = 1:length(DATA.helpfiles)
        uimenu(hm,'Label',DATA.helpfiles(j).label,'Callback',{@ShowHelp, DATA.helpfiles(j).filename});
    end
   uimenu(hm,'Label',sprintf('Version %s',strrep(DATA.vergversion,'verg.','')));
   if isfield(DATA.binoc{1},'ve') && ~isempty(DATA.binoc{1}.ve)
       v = strrep(DATA.binoc{1}.ve,'binoclean.','');
       uimenu(hm,'Label',sprintf('Binoc Version %s',v))
   end
 
  function OptionMenu(a,b,tag)     
      DATA = GetDataFromFig(a);
      on = get(a,'checked');
      f = get(a,'tag')
      if strcmp(on,'on');
          set(a,'checked','off');
          DATA.optionflags.(f) = 0;
      else
          set(a,'checked','on');
          DATA.optionflags.(f) = 1;
      end
      SendCode(DATA, 'optionflag');
      SetGui(DATA);
   
function ReBuildQuickMenu(DATA)
    hm = findobj(allchild(DATA.toplevel),'flat','Label','Quick','Tag','QuickMenu');
    delete(allchild(hm));
    BuildQuickMenu(DATA, hm);
          
  function BuildQuickMenu(DATA, hm)
        
    if isfield(DATA.quickexpts,'submenu')
        sms = {DATA.quickexpts.submenu};
        
%add submenus in the order they were created in the file        
    subs = unique(sms);
    if length(subs) > 1
        for j= 2:length(subs)
            id = find(strcmp(subs{j},sms));
            order(j) = id(1);
        end
        [a,b] = sort(order);
        if isempty(subs{b(1)})
            subs = subs(b(2:end)); %1st element is blank
        end
    end
    for j = 1:length(subs)
        sm(j) = uimenu(hm,'Label',subs{j});
    end
    for j = 1:length(DATA.quickexpts)
        k = find(strcmp(DATA.quickexpts(j).submenu,subs));
        if ~isempty(k)
            uimenu(sm(k),'Label',DATA.quickexpts(j).name,'Callback',{@verg, 'quick', DATA.quickexpts(j).filename});
        else
            uimenu(hm,'Label',DATA.quickexpts(j).name,'Callback',{@verg, 'quick', DATA.quickexpts(j).filename});
        end
    end
    end
    uimenu(hm,'Label','Today','Tag','TodayQuickList');
    uimenu(hm,'Label','Add...','Tag','AddToQuickList','callback',{@AddQuickMenu});

function AddQuickMenu(a,b)
    DATA = GetDataFromFig(a);
    [name, path] = uigetfile(['/local/' DATA.binoc{1}.monkey '/stims/*.*']);
    j = length(DATA.quickexpts)+1;
    DATA.quickexpts(j).name = name;
    DATA.quickexpts(j).filename = [path '/' name];
    if name
        hm = get(a,'parent');
        uimenu(hm,'Label',DATA.quickexpts(j).name,'Callback',{@verg, 'quick', DATA.quickexpts(j).filename});
    end
    
function SetElectrode(a,b, ei)
    DATA = GetDataFromFig(a);
    DATA.electrodeid = ei;
    SetMenuCheck(a,[], ei);
    SendCode(DATA, 'Electrode');
    set(DATA.toplevel,'UserData',DATA);
    
function DATA = RestartBinoc(DATA)
    if DATA.newbinoc == 1 %Don't restart if its new
        return;
    end
    fprintf('Restarting Binoc at %s\n',datestr(now));
    outprintf(DATA,'\\quit\n');
    if DATA.autoreopen == 0
        DATA.autoreopen = 2;
    end
    SetData(DATA);
    a = 1;
    tries = 0;
    while a ~= 0 && tries < 5
        [a,b] = system('open /local/bin/binoclean.app');
        if a ~= 0
            fprintf('Binoc wont open. Trying again in 1 sec\n');
        end
        pause(1);
        tries = tries+1;
    end
    
        
function MenuHit(a,b, arg)
    DATA = GetDataFromFig(a);
    if strcmp(arg,'bothclose')
        outprintf(DATA,'\\quit\n');
        [ok, reason] = CheckDayEnd(DATA);
        if ok
            CopyLog(DATA,'online');
            CopyLog(DATA,'penlog');
            CopyLog(DATA, 'bnc');
        elseif ~sum(strncmp(reason,{'Human'},6))
            fprintf('Not Copying Recording Logs beacuse %s\n',reason);
        end
        ExitVerg(DATA);
    elseif strcmp(arg,'restartbinoc')
        RestartBinoc(DATA);
    elseif strcmp(arg,'copylogs')
        ok = CheckDayEnd(DATA);
        if DATA.optionflags.py
            CopyLog(DATA, 'psych');
            CopyLog(DATA, 'serialpsych');
        else
            CopyLog(DATA, 'online');
            CopyLog(DATA, 'penlog');
            CopyLog(DATA, 'bnc');
        end
    elseif strncmp(arg,'checkdur',8)
        CheckTrialDurations(DATA,arg(9:end));
    elseif strcmp(arg,'choosefont')
        fn = uisetfont;
        DATA.font = fn;
        set(DATA.toplevel,'UserData',DATA);
    elseif strcmp(arg,'freereward')
        outprintf(DATA,'freerwd\n');
    elseif strcmp(arg,'updatepsych')
        UpdatePsychFile(DATA.binoc{1}.psychfile);
    elseif strcmp(arg,'onetrial')
        outprintf(DATA,'!onetrial\n');
    elseif strcmp(arg,'pipelog')       
        system([GetFilePath('perl') '/pipelog end']); %first stop any existing
        [a, prefix] = fileparts(DATA.binoc{1}.psychfile);
        prefix = regexprep(prefix,'[0-9][0-9][A-Z][a-z][a-z]20[0-9][0-9]','');
        prefix = regexprep(prefix,'DATE$','');
        srcfile = DATA.binoc{1}.psychfile;
        if strcmp(srcfile,'NotSet')
            srcfile = GetValue(DATA,'daylog');
            if ~exist(srcfile)
                vergwarning(sprintf('Cannot Start pipelog until Log exists. Come back after running a trial',srcfile));
            end
            system([GetFilePath('perl') '/pipelog ' DATA.binoc{1}.monkey ' log &']);
        else
            if ~exist(srcfile)
                vergwarning(sprintf('Cannot Start pipelog until %s exists. Come back after running a trial',srcfile));
            end
            system([GetFilePath('perl') '/pipelog ' DATA.binoc{1}.monkey ' ' prefix ' &']);
        end
        DATA.pipelog = 1;
        DATA = AddTextToGui(DATA,['Pipelog ' prefix]);
    elseif strcmp(arg,'setshake')
        outprintf(DATA,'usershake\n');
        
    elseif strcmp(arg,'showpsych')
        DATA = SetFigure('VergPsych', DATA);
        DATA.psych.blockmode = 'Current';
        DATA = PlotPsych(DATA, 'force');
        if ~isempty(DATA.err)
            fprintf('%s\n',DATA.err);
        end
        set(DATA.toplevel,'UserData',DATA);
    end
    
    
function CopyLog(DATA,type)
        
    
    dfile = strrep(DATA.binoc{1}.uf,'\','/');
    dfile = regexprep(dfile,'^[A-Z]:','');
    [a,b] = fileparts(dfile);
    fname=b;
    [a, fdir] = fileparts(a);
        
    if strcmp(type,'bnc')
        if strncmp('/local',DATA.binoc{1}.netpref,5) %only copy bnc file if its local
            bncfile = ['/local/' dfile '.bnc'];
            d = dir(bncfile);
            if length(d) == 1 && now - d.datenum < 1
                [a,b,c,x] = GetMonkeyName(bncfile);
                tgt = ['/b/data/'  x];
                if exist(tgt)
                    a = dir(tgt);
                    msg = sprintf('%s already exists: size %d\n %s size is %d\nCopy %s to %s?',tgt,a.size,bncfile,d.size.bncfile,tgt)
                else
                    msg = sprintf('Copy %s to %s?',bncfile,tgt)
                end
                if confirm(msg);
                    fprintf('Copying %s to %s\n',bncfile,tgt);
                    try
                        copyfile(bncfile,tgt);
                    end
                end
            end
        end
    elseif strcmp(type,'psych')
        src = DATA.binoc{1}.psychfile;
        [a,b] = fileparts(src);
        tgt = [DATA.binoc{1}.netpref '/' DATA.binoc{1}.monkey];
        if exist(tgt)
            tgt = [tgt '/' b];
            BackupFile(tgt,'print');
            fprintf('Copying %s to %s\n',src,tgt);
                try
                    copyfile(src,tgt);
                end
        else
            msg = sprintf('No Directory %s',tgt);
            acknowledge(msg);
        end
    elseif strcmp(type,'serialpsych')
        logfile = [DATA.cwd '/' DATA.binoc{1}.uf];
        tgt = [DATA.binoc{1}.netpref '/' DATA.binoc{1}.monkey '/serial/' DATA.binoc{1}.uf];
        fprintf('Copying %s to %s\n',logfile,tgt);
        try
            copyfile(logfile,tgt);
        end
    elseif strcmp(type,'online')
        logfile = ['/local/' DATA.binoc{1}.monkey '/' fname];
        d = dir(logfile);
        if length(d) == 1 && now - d.datenum < 1
            tgt = [DATA.binoc{1}.netpref '/' fname '.online'];
            ntgt = sprintf('/b/data/%s/%s/%s.online',DATA.binoc{1}.monkey, fdir, fname)
            if exist(tgt)
                msg = sprintf('Overwrite %s with %s?',tgt,logfile);           %s/d     
            else
                msg = sprintf('Copy %s to %s?(Network)\nor%s(PC)',logfile,ntgt,tgt);
            end
            
            T = 'Copying Log to Nework';
            ok = 0;
            strs = {'To Network' 'To PC' 'I''ll copy manually'};
            yn = myquestdlg(['Copy ' logfile ' to?'] , T, strs{:}, strs{1},DATA.font);
            ok = find(strcmp(yn,strs));
            if ok == 1
                tgt = ntgt;
            end
            if ok ~= 3
                fprintf('Copying %s to %s\n',logfile,tgt);
                try
                    copyfile(logfile,tgt);
                end
            end
        end
    elseif strcmp(type,'penlog')
        logfile = sprintf('/local/%s/pen%d.log', DATA.binoc{1}.monkey, DATA.binoc{1}.Pn(1));
        d = dir(logfile);
        if length(d) == 1 && now - d.datenum < 1
            tgtdir = sprintf('/b/data/%s/pens',DATA.binoc{1}.monkey);
            if ~exist(tgtdir)
                tgtdir = sprintf('/b/bgc/anal/%s',DATA.binoc{1}.monkey);
            end
            tgt = sprintf('%s/pen%d.log', tgtdir, DATA.binoc{1}.Pn(1));
            if exist(tgt)
                go = confirm(sprintf('%s Already Exists. Overwrite with  %s?',tgt,logfile),DATA.font);
            else
                go = confirm(sprintf('Copy %s to %s?',logfile,tgt),DATA.font);
            end
            if go
                fprintf('Copying %s to %s\n',logfile,tgt);
                try
                    copyfile(logfile,tgt);
                end
            end
        end
    end

 function SetExpt(a,b, type)
     DATA = GetDataFromFig(a);
     val = get(a,'value');
     str = get(a,'string');
    if strcmp(type,'et')
        outprintf(DATA,'et=%s\n',DATA.comcodes(DATA.expmenuvals{1}(val)).code);
    elseif strcmp(type,'e2')
        outprintf(DATA,'e2=%s\n',DATA.comcodes(DATA.expmenuvals{2}(val)).code);
    elseif strcmp(type,'e3')
        outprintf(DATA,'e3=%s\n',DATA.comcodes(DATA.expmenuvals{3}(val)).code);
    elseif strcmp(type,'fsd')
        if iscellstr(str)
        outprintf(DATA,'xyfsd=%s\n',str{val});
    else
        outprintf(DATA,'xyfsd=%s\n',str(val,:));
        end
        DATA.binoc{1}.xyfsd = val;
    elseif strmatch(type,{'st' 'bs'})
        id = strmatch(type,{'st' 'bs'});
        if val > 0
        DATA.stimtype(id) = val;
        outprintf(DATA,'mo=%s\n',DATA.stimtypenames{id});
        outprintf(DATA,'st=%s\n',DATA.stimulusnames{val});
        DATA.currentstim = id;
        SetGui(DATA);
        else
            fprintf('Cant set stim < 1\n');
        end
    end
    set(DATA.toplevel,'UserData',DATA);

    
function DATA = ReloadInitialSettings(DATA, name, varargin)
%DATA = ReloadInitialSettings(DATA, varargin)
% reads state variables from a stim file, except for 
% things that are day/session specicic

 xs = DATA.Header.dailyvars;

  strs = scanlines(name);
  
  for j = 1:length(xs)
      id = find(~strncmp(xs{j},strs,length(xs{j})));
      strs = strs(id);
  end
  DATA = ReadExptLines(DATA,strs,'file');
    
function DATA = LoadLastSettings(DATA, varargin)
%DATA = LoadLastSettings(DATA, varargin)
% reads state variables from expt backup files to restosr
% state - including se, id, penetraion details, ed.
% Only does this if the backup is less than half an hour old
% ...,'force')  overrides time test
    interactive = 0;
    j = 1;
    go  = 0;
    while j <= length(varargin)
        if strncmpi(varargin{j},'force',5)
            go = 1;
            yn = 'Yes';
        elseif strncmpi(varargin{j},'interactive',5)
            interactive = 1;
        end
        j = j+1;
    end
    
    if DATA.frombinocfid > 0 || DATA.verbose(4)
    end
        rfile = ['/local/' DATA.binoc{1}.monkey '/lean*.stm'];
        d = mydir(rfile);
        if DATA.frombinocfid > 0 || DATA.verbose(4)
            fprintf('%s has %d current expt files\n',rfile,length(d));
        end
        if isempty(d)
            return;
        end
        [a,id] = max([d.datenum]);
        d = d(id);
        if now - d.datenum < 0.5/24 %less than an half an hour old - read in settings
            go = 1;
            if interactive 
                cprintf('red', 'Looks Like restart.  Answer dialog\n');
                yn = questdlg(sprintf('Looks like you are re-starting Mid expt. Do you want to reload id/se from %s',d.filename),'Binoc has been working','Reload minimum','Reload all','No','Reload minimum');
                if ~strncmp(yn,'Reload',5)
                    go = 0;
                end
            end
        end
        if go
            txt = scanlines(d.name);
            if strcmp(yn,'Reload all')
                for j = 1:length(txt)
                    DATA = InterpretLine(DATA, txt{j},'tobinoc');
                end
            else
                for s = DATA.Header.dailyvars
                    id = find(strncmp(s,txt,length(s{1})));
                if ~isempty(id)
                    cprintf('blue','Setting %s from %s\n',txt{id(1)},d.name);
                    DATA = InterpretLine(DATA, txt{id(1)},'tobinoc');
                end
                end
            end
            if GetValue(DATA,'Pn') > 0
                outprintf(DATA,'!openpen\n');
            end
            DATA = GetState(DATA);
        end
    
function RecoverFile(a, b, type)
    DATA = GetDataFromFig(a);
%        fprintf('Recover called with %s\n',type);
    if strmatch(type,{'list'});
        rfile = ['/local/' DATA.binoc{1}.monkey '/lean*.stm'];
        d = dir(rfile);

        if isempty(d)
            return;
        end
        if strcmp(type,'list')
        hm = get(a,'parent');
        else
            hm = a;
        end
        c = get(hm,'Children');
        delete(c);
        uimenu(hm,'Label','List','callback',{@RecoverFile, 'list'});
        [a,id] = sort([d.datenum]);
        d = d(id);
        uimenu(hm,'Label',['Just Restore Pen/Expt Settings from ' d(end).name],'callback',{@RecoverFile, 'loadlast'});
        for j = 1:length(d)
            uimenu(hm,'Label',[d(j).name d(j).date(12:end)],'callback',{@RecoverFile, d(j).name(5:end)});
        end
        if isfield(DATA.state,'stimfile')
            uimenu(hm,'Label',['Restart and Reload from ' DATA.state.stimfile],'callback',{@RecoverFile, 'reloadstart'});
        end
        
    elseif strmatch(type,{'eo.stm' 'eb.stm' '0.stm' '1.stm' '2.stm' '3.stm' '4.stm' '5.stm'})
        rfile = ['/local/' DATA.binoc{1}.monkey '/lean' type];
        dfile = ['/local/' DATA.binoc{1}.monkey '/lean.today'];
        copyfile(rfile,dfile);
        DATA = ReadStimFile(DATA, dfile);
        set(DATA.toplevel,'UserData',DATA);
    elseif strcmp(type,'reloadstart')
        DATA = RestartBinoc(DATA);
        DATA = ReloadInitialSettings(DATA, DATA.state.stimfile);
    elseif strcmp(type,'loadlast')
        DATA = LoadLastSettings(DATA,'noninteractive','force');
        set(DATA.toplevel,'UserData',DATA);
    else
        if DATA.verbose(4)
            fprintf('Recover called with %s\n',type);
        end
    end
    
    
function ReadSlot(a,b,id)
    DATA = GetDataFromFig(a);
    quickname = sprintf('/local/%s/q%dexp.stm', DATA.binoc{1}.monkey,id);
    DATA = ReadStimFile(DATA, quickname);
    set(DATA.toplevel,'UserData',DATA);
    

    
function SaveSlot(a,b,id)
    DATA = GetDataFromFig(a);
    lb = sprintf('Expt%d: %s %s',id,DATA.stimulusnames{DATA.stimtype(1)},DATA.exptype{1});
    if ~strcmp(DATA.exptype{2},'e0')
        lb = [lb 'X' DATA.exptype{2}];
    end
    if ~strcmp(DATA.exptype{3},'e0')
        lb = [lb 'X' DATA.exptype{3}];
    end
    if DATA.optionflags.fS
        lb = [lb 'RC'];
    end
    set(a,'Label',lb);
    outprintf(DATA,'quicksave%d\n',id);
    AddTodayMenu(DATA, id, lb);

function AddTodayMenu(DATA, id,label)
    
    tag = sprintf('TodaySlot%d',id);
    it = findobj(DATA.toplevel,'Tag',tag);
    if isempty(it)
        it = findobj(DATA.toplevel,'Tag','TodayQuickList');
        if length(it) == 1
            uimenu(it,'Label',label,'callback',{@ReadSlot, id},'Tag',tag);
        end
    else
        set(it,'Label',label);
    end
    
function PushChanges(DATA)
    vpath = fileparts(mfilename('fullpath'));
    tgt = DATA.netmatdir;
    CheckFileUpdate([vpath '/verg.m'],[tgt '/verg.m'],'nobackup',DATA);
    CheckFileUpdate([vpath '/ServoDrive.m'],[tgt '/ServoDrive.m'],'nobackup',DATA);

function CheckForUpdate(DATA)
    
    vpath = mfilename('fullpath');
    tgt = DATA.localmatdir;
    if ~strncmp(vpath,DATA.localmatdir,length(DATA.localmatdir)) && ~strncmp(vpath,'/local/c/binoclean/matlab',25)
        str = sprintf('verg source is in %s, not local source directory (%s) Update to Current source directory?\n',fileparts(vpath),DATA.localmatdir);
        yn = myquestdlg(str,'Update Check','Yes','No','','Yes',DATA.font);
        if strcmp(yn,'Yes')
            tgt = fileparts(vpath);
        end
    end
    CheckFileUpdate([DATA.netmatdir '/verg.m'],[tgt '/verg.m'],'backup',DATA);
    CheckFileUpdate([DATA.netmatdir '/helpstrings.txt'],[tgt '/helpstrings.txt'],'change',DATA);
    CheckFileUpdate([DATA.netmatdir '/DownArrow.mat'],[tgt '/DownArrow.mat'],'new',DATA);
    CheckFileUpdate([DATA.netmatdir '/vergversion.m'],[tgt '/vergversion.m'],'change',DATA);
    CheckFileUpdate([DATA.netmatdir '/ServoDrive.m'],[tgt '/ServoDrive.m'],'backup',DATA);
    
 function CheckFileUpdate(src, tgt, chkmode,DATA)    
     if nargin < 3
        chkmode = 'change';
        DATA.font = [];
     end
    a = dir(src);
    b = dir(tgt);
    if isempty(a)
        fprintf('Missing %s - cant check update\n',src);
        return;
    end
    if isempty(b) %target does not exist - just copy
            try  %This often produces permission error
                copyfile(src,tgt);
            catch ME
                cprintf('errors',ME.message);
                fprintf('Error copying %s\n',tgt);
            end
    end
    if strncmp(chkmode,'new',3)         
        return;
    end
    if ~isempty(a) && ~isempty(b) && a.datenum > b.datenum
        yn = myquestdlg(sprintf('%s is newer. Copy to %s?',src,tgt),'Update Check','Yes','Compare','No','Yes',DATA.font);
        if strcmp(yn,'Yes')
            try  %This will produce and error becuase verg.m is in use. But the copy succeeds
                if strncmp(chkmode,'backup',3)
                    BackupFile(tgt);
                end

                [a,b,c] = copyfile(src,tgt);
            catch ME
                cprintf('errors',ME.message);
                fprintf('possible error copying %s\n',tgt);
            end
        elseif strcmp(yn,'Compare')
            visdiff(src,tgt);
        end
    end
    
    
function SaveFile(a,b,type)

    DATA = GetDataFromFig(a);
    if ~isfield(DATA,'stimfilename')
        DATA.stimfilename = ['/local/' DATA.binoc{1}.monkey '/stims/auto.stm'];
    end
    if strcmp(type,'current')
        filename = DATA.stimfilename;
        if exist(filename,'file')
            bakfile = [filename '.bak'];
            a = copyfile(filename,bakfile)
            if a
                fprintf('Saved backup in %s\n',bakfile);
            else
                fprintf('Couldnt save backup %s\n',bakfile)
            end
        end
        SaveExpt(DATA, filename);
    elseif strcmp(type,'fix')
        CheckStimFile(DATA,'new');
    elseif strcmp(type,'push')
        PushChanges(DATA)
    elseif strcmp(type,'update')
        CheckForUpdate(DATA)
    elseif strcmp(type,'layout')
        SaveLayout(DATA, DATA.layoutfile);
    elseif strcmp(type,'saveas')
        [a,b,c] = fileparts(DATA.stimfilename);
        if isempty(c)
            DATA.stimfilename = [DATA.stimfilename '.stm'];
        end
        [outname, path] = uiputfile(DATA.stimfilename, 'Save Expt As');
        if outname ~= 0
        SaveExpt(DATA, [path '/' outname]);
        end
    end
    
function SaveLayout(DATA, name)  

    fid = fopen(name,'w');
    for j = 1:length(DATA.winpos)
        it = findobj('tag',DATA.windownames{j},'type','figure');
        if length(it) == 1
            DATA.winpos{j} = get(it,'position');
        end
        fprintf(fid,'%s=%s\n',DATA.windownames{j},sprintf('%d ',DATA.winpos{j}));
    end
    fclose(fid);

function CheckTrialDurations(DATA, varargin)
%Never Essential - don't let this kill timer
    try
    T = DATA.Trials;
    plottype = 'none';
    j = 1;
    while j <= length(varargin)
        if strcmp(varargin{j},'EXPTOVER') && ~isempty(DATA.Expts)
            tid = DATA.Expts{end}.first:DATA.Expts{end}.last;
            T = DATA.Trials(tid);
        elseif sum(strcmp(varargin{j},{'hist' 'seq' 'expt'}))
            plottype = varargin{j};
        end
        j = j+1;
    end
    
    if isfield(T,'Nf') && isfield(T,'dur')
        Nf = [T.Nf];
        nf = [T.nf];
        durs = [T.dur];
        id = find(Nf == nf+1); %completed.
        err = 1+ durs(id).*DATA.binoc{1}.fz - Nf(id); %number of extra video frames
        if strcmp(plottype,'expt')
            xid = find([T(id).Start] > DATA.Expt.Start);
            err = err(xid);
        end
        if length(err) > 10
            x = sort(err);
            maxerr = round(x(end-3)); %allow a few wild ones
            if maxerr < 10
                maxerr = 10;
            end
        else 
            maxerr = 10;
        end
        if ~strcmp(plottype,'none')
            GetFigure(DATA.tag.plotwin,'parent',DATA.toplevels);
            hold off;
        end
        if strcmp(plottype,'expt')
            hist(err,[-10:0.5:maxerr]);
            title(sprintf('Real-Expected Duration at %.1fHz since %s',DATA.binoc{1}.fz,datestr(DATA.Expt.Start)));
            xlabel('frames');
        elseif strcmp(plottype,'hist')
            hist(err,[-10:0.5:10]);
            title(sprintf('Real-Expected Duration at %.1fHz',DATA.binoc{1}.fz));
            xlabel('frames');
        elseif strcmp(plottype,'seq')
           nomdur = Nf(id)./DATA.binoc{1}.fz;
           plot([T(id).Start],nomdur,'-');
           hold on;
           plot([T(id).Start],durs(id),'ro');
           datetick('x','HH:MM');
           xlabel('Time');                    
           ylabel('Duration (Sec)');
           title(sprintf('%d completed stim',length(id)));
        end
        if sum(err > 1) > length(err)/10 %10% trials are bad
            vergwarning(sprintf('%d/%d trials were too long',sum(err > 1),length(err)),'tellbinoc');
            eid = find(err > 1);
            for j = 1:length(eid)
                iT = T(id(eid(j)));
                fprintf('%s id %d:nf%d %.3f\n',datestr(iT.Start),iT.id,iT.nf,iT.dur);
            end
        end
        fprintf('%d/%d trials were too long\n',sum(err > 1),length(err));
    end
catch ME
CheckExceptions(ME);
end
    
function TextCallback(a,b)
    SendManualVals(a);

    
function DATA = SendManualVals(a, b)
% a can be a handel to the list, or
% a is DATA and b is the tag

 DATA = GetDataFromFig(a);

   if ishandle(a)
        tag = get(a,'Tag');
   elseif nargin > 1
       tag = b;
   else
       tag = 'All'; %default is to check all
   end
   if sum(strcmp(tag,{'All' 'ManualValuesButton'}))
       if DATA.listmodified(1)
           SendManualVals(DATA,'Expt1StimList');
       end
       if DATA.listmodified(2)
           SendManualVals(DATA,'Expt2StimList');
       end
       if DATA.listmodified(3)
           SendManualVals(DATA,'Expt3StimList');
       end
       DATA.listmodified = [0 0 0];
       SetData(DATA);
       if strcmp(tag, 'ManualValuesButton');
           delete(a);
       end
       return;
   end
   
    if ~ishandle(a) &&  nargin > 1
        a = findobj(allchild(DATA.toplevel),'flat','Tag',tag,'type','uicontrol');
    end
    str = get(a,'string');
    if strcmp(tag,'Expt2StimList')
        c = 'EB';
        o = 1+DATA.nextras(2);
        if size(str,1) < DATA.binoc{1}.n2+DATA.nextras(2)
            outprintf(DATA,'n2=%d\n',size(str,1)-DATA.nextras(2));
        end
    elseif strcmp(tag,'Expt3StimList')
        c = 'EC';
        o = 1+DATA.nextras(3);
    else
        o = 1+DATA.nextras(1);
        checknr =0;
        if size(str,1) > DATA.binoc{1}.nt+DATA.nextras(1) && checknr
            yn = questdlg('Nstim mismatch','That will Change N stims','OK','Cancel','OK');
        end
        if size(str,1) < DATA.binoc{1}.nt+DATA.nextras(1)
            outprintf(DATA,'nt=%d\n',size(str,1)-DATA.nextras(1));
        end
        c = 'EA';
    end

    for j = o:length(str)
            outprintf(DATA,'%s%d=%s\n',c,j-o,str{j});
    end
    outprintf(DATA,'EDONE\n');
     
    
function StimulusSelectMenu(a,b, fcn)
    DATA = GetDataFromFig(a);
    txt = ['mo=' fcn];
    DATA = InterpretLine(DATA, txt,'tobinoc');
    SetGui(DATA);
    set(DATA.toplevel,'UserData',DATA);
    
function EditValsMenu(a,b, fcn)
if strncmpi(fcn,'apply',5)
    SendManualVals(get(get(a,'parent'),'UserData'));
elseif strcmpi(fcn,'cancel')
    f = get(a,'parent');
    DATA = GetDataFromFig(get(f,'parent'));
    SendCode(DATA,'nt');
end
    
function EditText(a,b)
     sendvals = 0;
    DATA = get(gcf,'UserData');  %this only called from main window
    tag = get(a,'tag');
    e = find(strcmp(tag,{'Expt1StimList' 'Expt2StimList' 'Expt3StimList'}));
    if sum(strcmp(b.Key,{'uparrow' 'downarrow' 'leftarrow' 'rightarrow'})) == 0
        DATA.listmodified(e) = 1;
        h = getappdata(DATA.toplevel,'ApplyButton');
        if isempty(h) || ~ishandle(h)
            h = uicontrol(DATA.toplevel,'style','pushbutton','String','Apply (^G)','callback',@SendManualVals,'Tag','ManualValuesButton','units','normalized');
            setappdata(DATA.toplevel,'ApplyButton',h);
            x = get(h,'position');
            y = get(a,'position');
            x(1) = y(1);
            x(2) = y(2)+y(4);
            x(3) =1/5;
            set(h,'position',x);
        end
    end
     
     if double(b.Character) == 13
         str = get(a,'string');
     end
     if strcmp(b.Modifier,'control')
         if b.Key == 'g'
            sendvals = 1;
         end
     end
     if strcmp(b.Key,'escape')
         SendCode(DATA,'nt'); %makes binoc resend expt vals
     end
     
     if sendvals
         delete(h);
         rmappdata(DATA.toplevel,'ApplyButton');
%%MAJOR KLUDGE.  value of 'string' is NOT up do date inside the callback
% must force focus elsewere to update.  Poping up a window seems to work
         h = waitbar(0,'Applying Changes');
         close(h);
         SendManualVals(a);
         DATA = DrainBinocPipe(DATA);
%         SetExptItems(DATA); %should get done when receive 'EDONE' 
    end
     set(DATA.toplevel,'UserData',DATA);
    
   
 function EditList(a,b)
     id = get(a,'value')
     str = get(a,'string');
     c = get(a,'UserData');
     sendvals = 0;
%c(1) is index of current char
%c(2) is current line
%nb if this is changed with the mouse, watch out!
     DATA = GetDataFromFig(get(a,'parent'));
     
     
     
     if length(c) > 1 && c(2) ~= id  %line #changed
         c(1) = 0;
     end
     if isempty(c) || c(1) == 0;
         c= 1;
         str{id} = ' ';
     else
         c(1) = c(1)+1;
         if length(c) > 1
             set(a,'value',c(2));
             id = c(2);
         end
     end
     if double(b.Character) == 13
         newstr = input('Edit Values','Expt 1', length(str),str);
         c(1) = 0;
         set(a,'UserData',[c(1)]);
         sendvals = 1;
     elseif strcmp(b.Key,'downarrow')
         c(1) = 0;
         set(a,'UserData',[c(1) (id+1)]);
         sendvals = 1;
     elseif strcmp(b.Key,'uparrow')
         c(1) = 0;
         set(a,'UserData',[c(1) (id-1)]);
         sendvals = 1;
     elseif ~isempty(b.Character)
             if strcmp(b.Key,'backspace')
                 str{id} = ' ';
             else
                str{id}(c(1)) = b.Character;
             end
             set(a,'string',str);
             set(a,'UserData',[c(1) id],'value',id);
     end

if sendvals
    tag = get(a,'Tag');
    if strcmp(tag,'Expt2StimList')
        c = 'EB';
    elseif strcmp(tag,'Expt3StimList')
        c = 'EC';
    else
        c = 'EA';
    end

    for j = 1:length(str)
            outprintf(DATA,'%s%d=%s\n',c,j-1,str{j});
    end
    outprintf(DATA,'EDONE\n');
end

     fprintf('%d:%s %d (C%d)\n',id,str{id},get(a,'value'),c(1));
     
     
function DATA = SetNewPenetration(DATA)
    d = dir([DATA.cwd '/pen*.log']);
    for j = 1:length(d)
        pe(j) = sscanf(d(j).name,'pen%d');
    end
    [pes, id] = sort(pe,'descend');
    fprintf('Recent Penetrations:\n');
    %may not be 3 pends here
    for j = min([length(id) 3]):-1:1
        pendata{j} = ReadPen([DATA.cwd '/' d(id(j)).name],'noplot');
        fprintf('%d: %.1f,%.1f',pendata{j}.num,pendata{j}.pos);
        files = unique(pendata{j}.files);
        for k = 1:length(files)
            fprintf(' %s',files{k});
        end
        if isfield(pendata{j},'Electrode')
            fprintf(' Electrode %s',pendata{j}.Electrode);
            DATA.binoc{1}.Electrode = pendata{j}.Electrode;
        end
        if isfield(pendata{j},'ePr')
            fprintf(' Tube out %.1f mm',pendata{j}.ePr);
            DATA.binoc{1}.ePr = pendata{j}.ePr;
        end
        if isfield(pendata{j},'hemi')
            fprintf(' %sHemi',pendata{j}.hemi);
            DATA.binoc{1}.hemi = pendata{j}.hemi;
        end
        if isfield(pendata{j},'coarsemm')
            fprintf(' at %.1f mm',pendata{j}.coarsemm);
            DATA.binoc{1}.coarsemm = pendata{j}.coarsemm;
        end
        fprintf('\n');
    end
    DATA.binoc{1}.Pn = max(pe(pe < 2000))+1;
    SetGui(DATA);
    

    
function MenuBarGui(a,b)
     DATA = GetDataFromFig(a);
     str = get(a,'label');
     tag = get(a,'Tag');
     switch tag
         case 'NewPen'
             DATA = SetNewPenetration(DATA);
     end
     set(DATA.toplevel,'UserData',DATA);
    
function MenuGui(a,b)
     DATA = GetDataFromFig(a);
     strs = get(a,'string');
     val = get(a,'value');
     if iscellstr(strs)
     str = strs{val};
     else
     str = strs(val,:);
     end
     tag = get(a,'Tag');
     switch tag
         case 'ElectrodeType'
             DATA.binoc{1}.Electrode = str;
             DATA.electrodeid = val;
         case 'NewPen'
             DATA = SetNewPenetration(DATA);
         case 'VisualArea'
             DATA.binoc{1}.Vn = str;
             SendCode(DATA,'Vn');
             myprintf(DATA.penid,'VisualArea %s\n',str);
         case 'Monkey'
             mnk = DATA.binoc{1}.monkey;
             DATA.binoc{1}.monkey = str;
             NewMonkey(DATA, mnk);
             SendCode(DATA,'monkey');             
     end
     set(DATA.toplevel,'UserData',DATA);

     
function NewMonkey(DATA, mnk)
%mnk is old monkey    
    if isfield(DATA.binoc{1},'psychfile')
        if ~isempty(strfind(DATA.binoc{1}.psychfile,mnk))
            DATA.binoc{1}.psychfile = strrep(DATA.binoc{1}.psychfile,mnk,DATA.binoc{1}.monkey);
            fprintf('Setting Psych Log to %s\n',DATA.binoc{1}.psychfile);
            SendCode(DATA,'psychfile');
        end
    end

function DATA = AddComment(DATA, str, src)
    
    DATA = AddStatusLine(DATA,str,'comment');
    DATA.Comments(end+1).comment = str;
    DATA.Comments(end).date = now;
    DATA.Comments(end).src = src;
    if isfield(DATA.binoc{1},'ed')
        DATA.Comments(end).ed = DATA.binoc{1}.ed;
    end
    CommentPopup(DATA,[],'update');
    if nargout == 0 %if ask for return, caller will Call SetData
        SetData(DATA);
    end
     
 function TextGui(a,b, type)
     DATA = GetDataFromFig(a);
     str = get(a,'string');
     if nargin == 2
         type = get(a,'Tag');
     end
     switch type
         case 'RptExpts'
             DATA.rptexpts = str2num(str);
             if isempty(DATA.rptexpts)
                 DATA.rptexpts = 0;
             end
             set(DATA.toplevel,'UserData',DATA);
         case 'cm'
             outprintf(DATA,'cm=%s\n',str);
             DATA = AddComment(DATA,str,'user');
             LogCommand(DATA,sprintf('cm=%s',str)); %sets data
             set(a,'string','');
         case 'nt'
             outprintf(DATA,'nt=%d\n',str2num(str));
             ReadFromBinoc(DATA);
         case 'n2'
             DATA.nstim(2) = str2num(str);
             outprintf(DATA,'n2=%d\n',DATA.nstim(2));
             ReadFromBinoc(DATA);
         case 'n3'
             DATA.nstim(3) = str2num(str);
             outprintf(DATA,'n3=%d\n',DATA.nstim(3));
             ReadFromBinoc(DATA);
         case 'em'
             DATA.mean(1) = str2num(str);
             outprintf(DATA,'em=%.8f\n',DATA.mean(1));
             ReadFromBinoc(DATA);
         case 'm2'
             DATA.mean(2) = str2num(str);
             outprintf(DATA,'m2=%.8f\n',DATA.mean(2));
             ReadFromBinoc(DATA);
         case 'm3'
             DATA.mean(3) = str2num(str);
             outprintf(DATA,'m3=%.8f\n',DATA.mean(3));
             ReadFromBinoc(DATA);
         case {'ei' 'i2' 'i3'}
             outprintf(DATA,'%s=%s\n',type,str);
             ReadFromBinoc(DATA);
         case 'st'
             DATA.stimtype(1) = strmatch(str,DATA.stimulusnames);
             set(DATA.toplevel,'UserData',DATA);
         case 'bs'
             DATA.stimtype(2) = strmatch(str,DATA.stimulusnames);
             set(DATA.toplevel,'UserData',DATA);
         case 'uf'
             DATA = InterpretLine(DATA, ['uf=' str],'fromgui','tobinoc');
             set(DATA.toplevel,'UserData',DATA);
             SetGui(DATA);
         otherwise
             DATA.binoc{DATA.currentstim}.(type) = str2num(str);
             outprintf(DATA,'%s=%s\n',type,str);
             ReadFromBinoc(DATA);
             
             
     end
             
 function outprintf(DATA,varargin)
 %send to binoc.  ? reomve comments  like in expt read? 
 show = 0;
 waitforbinoc = 1;
 if strcmp(varargin{1},'-show')
     show = 1;
     varargin = varargin(2:end);
 end
 [a,b] = cellstrcmp('-nowait',varargin);
 if a
     varargin = varargin(setdiff(1:length(varargin),find(b)));
     waitforbinoc = 0;
 end
 
 if DATA.network
     if DATA.binocisup
     str = sprintf(varargin{:});
     strs = split(str,'\n');
     for j = 1:length(strs)
         if ~isempty(strs{j})
%replace leading# with '//' so that it gets past http             
             if strs{j}(1) == '#' 
                 strs{j} = ['//' strs{j}(2:end)];
             end
                 
             str = [DATA.ip strs{j}];
             if DATA.verbose(1) || show
                 fprintf('%s\n',str);
             end
             if(DATA.tobinocfid > 0)
                 fprintf(DATA.tobinocfid,'%s:%s\n',datestr(now),str);
             end
             if waitforbinoc
                 ts = now;
                 [bstr, status] = urlread(str,'Timeout',2);
                 if ~isempty(bstr)
                     fprintf('Binoc replied with %s\n',bstr);
                 elseif DATA.verbose(5)
                     fprintf('Binoc returned in %.3f\n',mytoc(ts));
                 end
             end
         end
     end
     end
 elseif DATA.outid
     fprintf(DATA.outid,varargin{:});
 end
     
 function myprintf(fid,varargin)
     show = 0;
     j = 1;
     while j <= length(varargin)
         if strncmpi(varargin{j},'-show',5)
             varargin = varargin(j+1:end);
             show = 1;
         elseif strncmpi(varargin{j},'-date',5)
             varargin = varargin(j+1:end);
             if fid > 0
                 fprintf(fid,'%s ',datestr(now));
             end
         end
         j = j+1;
     end
     if fid > 0
         fprintf(fid,varargin{:});
     end
     if show
         fprintf(varargin{:});
     end
     
 function SendStr(a,b, str)
     DATA = GetDataFromFig(a);
     outprintf(DATA,'%s\n',str);
     DATA = ReadFromBinoc(DATA);
     SetGui(DATA);
    
 function CheckTimerHit(a,b, flag)
     DATA = GetDataFromFig(a);
     CheckTimer(DATA);
     PauseRead(DATA,0); 
 
 function DATA = ReopenPipes(a,b)
     DATA = ReadIO(a,b, 6);
     
 function DATA = ReadIO(a,b, flag)
     DATA = GetDataFromFig(a);

    if strcmp(flag,'openlog')
        DATA = OpenBinocLog(DATA,'frombinoc');
        SetData(DATA);
    elseif flag == 2
         ts = now;
         DATA = GetState(DATA,'ReadIO',1);
         fprintf('Manual State took %.2f',mytoc(ts));
         set(DATA.toplevel,'UserData',DATA);
         SetGui(DATA);
         fprintf('  +GUI setting %.2f\n',mytoc(ts));
         PauseRead(DATA,0); %force off here
     elseif flag == 3
         stop(DATA.timerobj)
         outprintf(DATA,'NewMatlab\n');
        DATA = ReadFromBinoc(DATA,'reset');   
        SetGui(DATA);
        start(DATA.timerobj);
     elseif flag == 4
        stop(DATA.timerobj)
     elseif flag == 5
        DATA = ReadFromBinoc(DATA,'reset');   
        outprintf(DATA,'\neventcontinue\nEDONE\n');
        setappdata(DATA.toplevel,'WaitingForDlg',0);
        if ~strcmp(get(DATA.timerobj,'Running'),'on')
        start(DATA.timerobj);
        end
     elseif flag == 6  % Reopen pipes
        stop(DATA.timerobj);
        if DATA.inexpt
             DATA.optionflags.do = 0;
         end
         DATA = OpenPipes(DATA, 0);
         SendState(DATA,'all');
         if GetValue(DATA,'Pn') > 0
             outprintf(DATA,'!openpen\n');
         end
         DATA = GetState(DATA,'Reopen');
         SetGui(DATA);
         StartTimer(DATA);
         set(DATA.toplevel,'UserData',DATA);
        
     elseif flag == 7 %obsolete
         if DATA.verbose(1) > 0
             DATA.verbose(1) = 0;
             set(a,'Label','verbose pipes');
         else
             DATA.verbose = 2;
             set(a,'Label','Quiet pipes');
             DATA = OpenBinocLog(DATA,'frombinoc');
         end
         outprintf(DATA,'verbose=%d\n',DATA.verbose(2));
         set(DATA.toplevel,'UserData',DATA);
     elseif flag == 8
        DATA = ReadFromBinoc(DATA,'reset','verbose2');
     else
        DATA = ReadFromBinoc(DATA);   
        SetGui(DATA);
     end
        CheckTimer(DATA);

function SetMenuToggle(a,b,flag)        
     DATA = GetDataFromFig(a);
     if ~isfield(DATA,flag)
         DATA.(flag)= 1;
     else
         DATA.(flag) = ~DATA.(flag);
     end
     if DATA.(flag)
         set(a,'checked','on');
     else
         set(a,'checked','off');
     end
     SetData(DATA);
     
function SetVerbose(a,b, flag)
     DATA = GetDataFromFig(a);
%verbose(4) Prints each trial outcome to console
%verbose(2) is sent to binoc controls prints/NSlogs there
%verbose(7) makes verg print what it receives from binoo

     if flag == 0
         DATA.verbose = zeros(size(DATA.verbose));
         SetMenuCheck(a,'exclusive', 1);
     elseif strcmp(flag,'openlogto')
         DATA= OpenBinocLog(DATA,'tobinoc')
         set(a,'checked','on');
     elseif strcmp(flag,'openlogfrom')
         DATA= OpenBinocLog(DATA,'frombinoc')
         set(a,'checked','on');
     elseif strcmp(flag,'nopopup')
         DATA.nowarning = ~DATA.nowarning;
         if DATA.nowarning
             vergwarning('-nopopup');
         else
             vergwarning('-popup');
         end
           
         SetMenuCheck(a, DATA.nowarning);
     else
         DATA.verbose(flag) = ~DATA.verbose(flag);
         if DATA.verbose(flag)
             set(a,'checked','on');
         else
             set(a,'checked','off');
         end             
         if flag ==2
                    outprintf(DATA,'verbose=%d\n',DATA.verbose(2));
         end

     end
     set(DATA.toplevel,'UserData',DATA);
        
 function DATA = CheckExptMenus(DATA)
     new = 0;
     for j = 1:3
         it = findobj(allchild(DATA.toplevel),'flat','Tag',['Expt' num2str(j) 'List']);
         id = strmatch(DATA.exptype{j},DATA.expmenucodes{j},'exact');
         mstr = get(it,'String');
         if isempty(id)
             DATA.expmenucodes{j} = {DATA.expmenucodes{j}{:} DATA.exptype{j}};
             a = find(strcmp(DATA.exptype{j},{DATA.comcodes.code}));
             if isempty(a)
                 str = 'Unknown';
             else
                 str = DATA.comcodes(a).label;
             end
             DATA.expstrs{j}{end+1} = str;
             set(it,'string',DATA.expstrs{j});
             new = new+1;             
         elseif id > length(mstr)
             fprintf('GUI Expt %d list was reset\n',j)
             set(it,'string',DATA.expstrs{j});
         end
     end
if new
    set(DATA.toplevel,'UserData',DATA);
end

function StartTimer(DATA, onoff)
    if nargin == 1
        onoff = 1;
    end
    if isfield(DATA,'timerobj') & isvalid(DATA.timerobj)
        on = strcmp(get(DATA.timerobj,'Running'),'on');
        if onoff == 0&& on
            %              stop(DATA.timerobj);
        elseif onoff == 1 && on == 0
            start(DATA.timerobj);
        end
    end
         
function current_state = PauseRead(DATA, onoff)
%current_state = PauseRead(DATA, onoff)
%onoff 1 = pause is ON.
%returns the state BEFORE the call. It always sets the desired state
if isfield(DATA,'toplevel') && isfigure(DATA.toplevel)
    current_state = getappdata(DATA.toplevel,'PauseReading');
    if isempty(current_state)
        current_state = 0;
    end
    if nargin == 2
        if onoff == current_state && current_state == 0
            fprintf('Pause is already %d\n',onoff);
        end
        setappdata(DATA.toplevel,'PauseReading',onoff);
        if isfield(DATA,'timerobj') & isvalid(DATA.timerobj)
            on = strcmp(get(DATA.timerobj,'Running'),'on');
            if onoff && on
  %              stop(DATA.timerobj);
            elseif onoff == 0 && on == 0
  %              start(DATA.timerobj);
            end
        end
    else
        onoff = current_state;
    end
elseif isfield(DATA,'pausereading')
%? need this. Timer should not start before figure is up?    
    current_state = DATA.pausereading;
    if nargin == 2
        DATA.pausereading = onoff;
    end    
    current_state = 0;
end

function SetExptItems(DATA, varargin)
    SetTextItem(DATA.toplevel,'Expt1Nstim',DATA.binoc{1}.nt);
    SetTextItem(DATA.toplevel,'Expt2Nstim',DATA.nstim(2));
    SetTextItem(DATA.toplevel,'Expt3Nstim',DATA.nstim(3));
    SetTextItem(DATA.toplevel,'Expt1Incr',DATA.binoc{1}.ei);
    SetTextItem(DATA.toplevel,'Expt2Incr',DATA.binoc{1}.i2);
    SetTextItem(DATA.toplevel,'Expt3Incr',DATA.binoc{1}.i3);
    SetTextItem(DATA.toplevel,'Expt1Mean',DATA.mean(1));
    SetTextItem(DATA.toplevel,'Expt2Mean',DATA.mean(2));
    SetTextItem(DATA.toplevel,'Expt3Mean',DATA.mean(3));
    SetTextItem(DATA.toplevel,'DataFileName',DATA.binoc{1}.uf);
    SetTextItem(DATA.toplevel,'binoc.nr',DATA.binoc{1}.nr);


 function SetGui(DATA,varargin)
     if ~isfield(DATA,'toplevel')
         return;
     end
     j = 1;
     while j <= length(varargin)
         if strncmpi(varargin{j},'ifneed',3)
             if isfield(DATA,'guiset') && DATA.guiset > 0
                 return;
             end
         elseif strncmpi(varargin{j},'set',3)
             set(DATA.toplevel,'UserData',DATA);
         end
         j = j+1;
     end
     paused = PauseRead(DATA); %get current state
     if paused ==0
        PauseRead(DATA,1);
     end
    DATA= CheckExptMenus(DATA);
    SetExptItems(DATA);
    SetTextItem(DATA.toplevel,'RptExpts',DATA.rptexpts);
    id = strmatch(DATA.exptype{1},DATA.expmenucodes{1},'exact');
    SetMenuItem(DATA.toplevel, 'Expt1List', id);
    id = strmatch(DATA.exptype{2},DATA.expmenucodes{2},'exact');
    SetMenuItem(DATA.toplevel, 'Expt2List', id);
    id = strmatch(DATA.exptype{3},DATA.expmenucodes{3},'exact');
    SetMenuItem(DATA.toplevel, 'Expt3List', id);
    SetMenuItem(DATA.toplevel, 'ForegroundType', DATA.stimtype(1));
    SetMenuItem(DATA.toplevel, 'BackgroundType', DATA.stimtype(2));
    SetToggleItem(DATA.toplevel, 'XYR', DATA.showxy(1));
    SetToggleItem(DATA.toplevel, 'XYL', DATA.showxy(2));
    SetToggleItem(DATA.toplevel, 'XYB', DATA.showxy(3));
    it= findobj(allchild(DATA.toplevel),'flat','Tag','CurrentStimLabel');
    set(it,'string',DATA.stimlabels{DATA.currentstim});
    if DATA.currentstim > 1
        set(it,'backgroundcolor','r');
    else
        set(it,'backgroundcolor', DATA.windowcolor);
    end
    it = findobj(allchild(DATA.toplevel),'flat','Tag','RunButton');
    if DATA.inexpt
        set(it,'string','Cancel');
    elseif DATA.optionflags.ts
        set(it,'string','Store');
        set(it,'backgroundcolor', DATA.windowcolor);
    else
        set(it,'string','Run','backgroundcolor',DATA.windowcolor);
        if isfield(DATA.optionflags,'py') && ~DATA.optionflags.py  %%Not Human Psych
            if strncmp(DATA.smrfile,DATA.binoc{1}.uf,length(DATA.binoc{1}.uf))
                set(it,'backgroundcolor','r');
            else
                set(it,'backgroundcolor', DATA.windowcolor);
            end
        end
    end
    
    it= findobj(allchild(DATA.toplevel),'flat','Tag','Monkey');
    if ~isempty(it)
        strs = get(it,'string');
        id = find(strcmp(DATA.binoc{1}.monkey,strs));
        if ~isempty(id)
            set(it,'value',id);
        end
    end
 
    
    ot = findobj(allchild(0),'flat','tag',DATA.windownames{2},'type','figure'); %options window
       f = fields(DATA.optionflags);
       otchild = allchild(ot);
       topchild = allchild(DATA.toplevel);
       for j = 1:length(f)
           if length(ot) == 1 %set value in options window
               it = findobj(otchild,'flat','Tag',f{j});
               set(it,'value',DATA.optionflags.(f{j}));
           end
           it = findobj(topchild,'flat','Tag',f{j},'type','uicontrol');
           if length(it) == 1
               set(it,'value',DATA.optionflags.(f{j}));
           end
       end
    ot = findobj(allchild(0),'flat','tag',DATA.windownames{3},'type','figure'); %software Offset
    SetSoftOffWindow(DATA,ot);
    
    
    ot = findobj(allchild(0),'flat','tag',DATA.windownames{9},'type','figure');
    if ~isempty(ot)
        otchild = allchild(ot);
        f = {'Pn' 'ePr' 'coarsemm'};
        for j = 1:length(f)
            it = findobj(otchild,'flat', 'tag',f{j},'type','uicontrol','style','edit');
            if length(it) ==1
                set(it,'string',sprintf('%d',GetValue(DATA,f{j})));
            end
        end
        it = findobj(otchild, 'flat','tag','ElectrodeType','type','uicontrol');
        if length(it) ==1
            estr = GetValue(DATA,'Electrode');
            strs = get(it,'string');
            sval = find(strcmp(estr,strs));
            if length(sval) ==1
                set(it,'value',sval);
            end
        end
        it = findobj(otchild,'flat', 'tag','hemisphere','type','uicontrol');
        if length(it) ==1
            estr = GetValue(DATA,'hemi');
            strs = cellstr(get(it,'string'));
            sval = find(strcmp(estr,strs));
            if length(sval) ==1
                set(it,'value',sval);
            end
        end
    end
    [a,j] = min(abs(DATA.binoc{1}.xyfsd - DATA.xyfsdvals));
    ot = findobj(allchild(DATA.toplevel),'flat','tag','FSD','type','uicontrol');
    set(ot,'value',j);

    ot = findobj(allchild(DATA.toplevel),'flat','tag','OptionContextMenu');
    if ~isempty(ot);
        m = get(ot,'children');
        for k = 1:length(m)
        c = get(m(k),'children');
        for j = 1:length(c)
            tag = get(c(j),'tag');
            if isfield(DATA.optionflags,tag) && DATA.optionflags.(tag)
                set(c(j),'checked','on');
            else
                set(c(j),'checked','off');
            end
        end
        end
    end

    ot = findobj(allchild(DATA.toplevel),'flat','tag','UffButton');
    if strncmp(DATA.smrfile,DATA.binoc{1}.uf,length(DATA.binoc{1}.uf))
        set(ot,'BackGroundColor',DATA.windowcolor);
    else
        set(ot,'BackGroundColor','r');
    end
    
    if DATA.verbose(4)
        fprintf('SetGUI: Ex%d\n',DATA.inexpt);
    end
    DATA.guiset = 1;

    
    if paused == 0 %only release pause if we weren't paused at the start.
        PauseRead(DATA,0);
    end

    
    

 function SetTextItem(top, tag, value, varargin)
 
     it = findobj(allchild(top),'flat','Tag',tag);
 if ~isempty(it)
     if ischar(value)
     set(it,'string',value);
     else
     set(it,'string',num2str(value));
     end
 end

 function SetToggleItem(top, tag, value, varargin)
 it = findobj(allchild(top),'flat','Tag',tag);
 if ~isempty(it)
     set(it,'value',value);
 end
 
function SetMenuItem(top, tag, value, varargin)
if length(value) == 1
     it = findobj(allchild(top),'flat','Tag',tag);
     if ~isempty(it)
         str = get(it,'string');
         if value > size(str,1)
             fprintf('Value %d out of range for %s\n',value, tag);
         else
             set(it,'value',value);
         end
     end
elseif length(value) > 1
     it = findobj(allchild(top),'flat','Tag',tag);
     if ~isempty(it)
         set(it,'value',value(1));
     end
end


function CheckServo(a,b, fig, varargin);
    DATA = get(fig,'UserData');
    if isfigure(DATA.servofig)
        ServoDrive('readposition','quiet');
    end

function CheckInput(a,b, fig, varargin)
    persistent lastread;
    DATA = get(fig,'UserData');
    
    pauseread = getappdata(DATA.toplevel, 'PauseReading');
    if pauseread && ~isempty(lastread)
        if DATA.verbose(1)
            fprintf('Paused...');
        end
        ts = mytoc(lastread);
        if ts > DATA.pausetimeout && DATA.pausetimeout > 0
            setappdata(DATA.toplevel,'PauseReading',0);
            vergwarning(sprintf('Paused for %.2f seconds. Unlocking.',ts));
        end
        return;
    end
    lastread = now;
    if DATA.timerperiod > 0.1 && DATA.tobinocfid > 0
        a = datevec(now);
        fprintf(DATA.tobinocfid,'Ping at %2d:%2d.%.3f\n',a(4),a(5),a(6));
    end
    if DATA.network == 2 %don't ping binoc during stims
        ts = (lastread - DATA.lastreadtime)/(60*60*24);
        if ts > 5
            fprintf('Have not heard from binoc for %.1f sec',ts);
        end
        a = dir('/tmp/binocstimisup');
        b = dir('/tmp/binocstimisdone');
        if ~isempty(a) && ~isempty(b) && a.datenum > b.datenum %binoc is in stim            
 %           myprintf(DATA.frombinocfid,'Binoc In Stim at %s\n',datestr(now));
            return;
        elseif isempty(b) && ~isempty(a) %in stim and done file deleted
%            myprintf(DATA.frombinocfid,'Binoc In Stim(D) at %s\n',datestr(now));
            return;
        end
    end
    try
        ReadFromBinoc(DATA, 'auto');
    catch ME
        cprintf('red','Error in timer call\n');
        s = getReport(ME);
        fprintf(s);
    end
    if DATA.confused
        fprintf('Getting State from binoc\n');
        DATA = GetState(get(DATA.toplevel,'UserData'),'Confused');
        DATA.confused = 0;
        SetData(DATA);
        fprintf('inexpt is now %d\n',DATA.inexpt);
    end
    if DATA.canceltimer > 0 
        if isfield(DATA,'Expts') && ~isempty(DATA.Expts)
        if DATA.inexpt
            if(now-DATA.Expts{end}.Start > (DATA.canceltimer ./ (24 * 60* 60)));            
                DATA = RunButton(DATA,'Cancel',1);
            end
        else
            if(now-DATA.Expts{end}.End > (DATA.canceltimer ./ (24 * 60* 60)));            
                DATA = RunButton(DATA,'Run',1);
            end            
        end
        end
    end
    if DATA.verbose(1) > 1
    fprintf('Timer read over at %s\n',datestr(now));
    end
    
    
 function [DATA, fullstr] = DrainBinocPipe(DATA, varargin)

waitforstim = 0;
ok = 1;
     
     j = 1; 
     while j <= length(varargin)
         if strncmpi(varargin{j},'waitforstim',7)
             waitforstim = 1;
             ok = 0;
         end
         j = j+1;
     end
     fullstr = '';
     if isfield(DATA,'toplevel')
         pausestate = getappdata(DATA.toplevel,'PauseReading');
     else
         pausestate = 1; % do nothing
     end
     if pausestate == 0
         PauseRead(DATA,1);
     end
     newmagic = DATA.binoc{1}.magic+1;
     outprintf(DATA,'magic=%d\n',newmagic);
     outprintf(DATA,'magic=');
     ts = now;
     trialcount = DATA.nt;
     while (DATA.binoc{1}.magic ~= newmagic || ok == 0) && mytoc(ts) < DATA.draintimeout
         [DATA, str] = ReadHttp(DATA);
         if DATA.nt > trialcount
             ok = 1;
         end
         fullstr = [fullstr str];
         myprintf(DATA.frombinocfid,'Drain at %s\n',datestr(now));
     end
     DATA.readdur = mytoc(ts);
     if DATA.readdur >= DATA.draintimeout
         fprintf('Read timed out after %.2f sec\n',DATA.readdur);
     end
     if pausestate == 0
         PauseRead(DATA,0);
     end
    
    
function [DATA, str] = ReadHttp(DATA, varargin)
j = 1;
expecting= 0;
expecttime = 1;
silent = 1;
str = [];
persistent httpbusy;
persistent lastts;
srcchr = '';
     
if isempty(lastts)
    lastts = 0;
end
     
if isempty(httpbusy)
    httpbusy = 0;
end
    args = {};
    while j <= length(varargin)
         if strncmpi(varargin{j},'verbose',5)
             verbose(1) = 1;
         elseif strncmpi(varargin{j},'expecting',5)
            expecting = 1;
            if length(varargin) > j && isnumeric(varargin{j+1})
                silent = 1;
                j = j+1;
                expecttime = varargin{j};
            end
         elseif strncmpi(varargin{j},'auto',4)
             srcchr = 'T';
         elseif strncmpi(varargin{j},'paused',4)
             args = {args{:} varargin{j}};
         elseif strncmpi(varargin{j},'reset',5)
             httpbusy = 0;
         end
         j = j+1;
    end

    DATA.httpbusy = httpbusy;
    DATA.readdur = 0;
     if httpbusy
         return;
     end
     
     ts = now;
    [str, status] = urlread([DATA.ip 'whatsup']);
    if expecting && DATA.verbose(4) %verg state
        fprintf('Expect: %d bytes\n',length(str));
    end
    while expecting && strncmp(str,'SENDING000000',13)
        httpbusy = 1;
        [str, status] = urlread([DATA.ip 'whatsup']);
        if silent == 0
        fprintf('0 Bytes, but expect input. Trying again.\n');
        if strncmp(str,'SENDING000000',13)
            fprintf('Nothing after %.2f\n',mytoc(ts));
        else
            fprintf('Second Try got %s\n',str(1:13));
        end
        end
        if ~isnumeric(expecttime)
            expecttime = 1;
        end
        if mytoc(ts) > expecttime
            expecting = 0;
            myprintf(DATA.frombinocfid,'-show','Still No Response from Binoc. I give up\n');
        end
     end

    took = mytoc(ts);
    DATA.readdur = took;
    if DATA.verbose(5)
        fprintf('Binoc status%d:%s\n',status,str);
        fprintf('Read took %.3f,%.3f\n',took,mytoc(ts));
    elseif DATA.verbose(7) && ~strncmp(str,'SENDING000000',13)
        fprintf('Binoc status%d:%s\n',status,str);
    end
    if isempty(str)
        [DATA,lastts] = CheckForNewBinoc(DATA,lastts);
    else
        DATA = InterpretLine(DATA,str,['frombinoc' srcchr],args{:});
        DATA.lastreadtime = now;
    end
%    myprintf(DATA.frombinocfid,'ReadBinoc%s:  %s',datestr(now),str);
    if ~isfield(DATA,'trialcounts')
    end
    lastts = ts;
httpbusy = 0;
 
     
     
 function [DATA, str] = ReadFromBinoc(DATA, varargin)
     global rbusy;
     persistent lastts;
     persistent lasttnone;
     
     if isempty(lastts)
         lastts = 0;
         lasttnone = 0;
     end
     verbose = DATA.verbose;
     autocall = 0;
     expecting = 0;
     str = [];

if DATA.network
    if DATA.binocisup
        [DATA, str] = ReadHttp(DATA, varargin{:});
        SetData(DATA);
    end
    return;
end
     
     j = 1;
     while j <= length(varargin)
         if strncmpi(varargin{j},'verbose',5)
             verbose(1) = 1;
             if strncmpi(varargin{j},'verbose2',8)
                 verbose = ones(size(verbose));
             end
         elseif strncmpi(varargin{j},'auto',4)
             autocall = 1;
         elseif strncmpi(varargin{j},'expecting',5)
            expecting = 1;
         elseif strncmpi(varargin{j},'reset',5)
             rbusy = 0;
         elseif ischar(varargin{j})
             fprintf('%s',varargin{j});
         end
         j = j+1;
     end
     if rbusy > 0
         fprintf('ReadBinoc busy at %s since %s\n',datestr(now),datestr(rbusy));
         return;
     end
     ts = now;
     rbusy = ts;    
         
     if DATA.outid <= 1
         return;
     end
     if verbose(7)
     fprintf('%s:',datestr(now,'HH:MM:SS.FFF'))
     end
     myprintf(DATA.frombinocfid,'ReadBinoc %.3f\n',mytoc(DATA.starttime));
     outprintf(DATA,'whatsup\n');
     a = fread(DATA.inid,14);
     if verbose(7)
         fprintf('OK\n');
     end
     if expecting && strcmp(char(a'),'SENDING000000')
         fprintf('0 Bytes, but expect input. Trying again.\n');
         a = fread(DATA.inid,14);
     end
     rbusy = 0;
     if DATA.frombinocfid > 0 && ~isempty(a)
         fprintf(DATA.frombinocfid,'%sB %s',datestr(now),char(a'));
     end
     if strncmp(char(a'),'SENDINGstart1',12)
        a = fread(DATA.inid,14);
        nbytes = sscanf(char(a'),'SENDING%d');
        fprintf('Start String\n');
     elseif strncmp(char(a'),'SENDING',7)
         nbytes = sscanf(char(a'),'SENDING%d');
         lastts = ts;
     else
         s = char(a');
         if ts - lasttnone > 5 /(24 * 60 * 60)
             fprintf('No Bytes %s %s\n',s,datestr(ts));
             lasttnone = ts;
         end
         if length(s)
             id = strfind(s,'SENDING')
             if length(id)
                 fprintf('Found SENDING at char %d\n',id);
                 nbytes = sscanf(s(id(1)),'SENDING%d')
             else
                 a = s(end);
                 while char(a) ~= 'G' | strcmp(s(end-6:end),'SENDING') == 0
                     a = fread(DATA.inid,1);
                     s = [s char(a)];
                 end
                 a = fread(DATA.inid,7);
                 nbytes = sscanf(char(a'),'%d');
                 fprintf('Read %s\n',s);
             end
         else
             nbytes = 0;
             [DATA, lastts] = CheckForNewBinoc(DATA, lastts);
         end
     end
     if verbose(7)
     fprintf('Need %d bytes\n',nbytes);
     end
     if nbytes > 0
         a = fread(DATA.inid,nbytes);
         if verbose(7)
         fprintf('%s',char(a'));
         fprintf('Read %d bytes took %.2f\n',length(a),mytoc(ts));
         end
         myprintf(DATA.frombinocfid,'Read %d bytes at %s\n',nbytes,datestr(now));
         myprintf(DATA.frombinocfid,'%s',char(a'));
         DATA = InterpretLine(DATA,char(a'),'frombinoc');
         if isfield(DATA,'toplevel')
             set(DATA.toplevel,'UserData',DATA);
         end
         myprintf(DATA.frombinocfid,'Finished at %s\n',datestr(now));
         lastts = ts;
     end
     rbusy = 0;

     
     
     
 function [DATA, lastts] = CheckForNewBinoc(DATA, lastts)

     go = 0;
if DATA.ready == 0
    return;
end
if nargin > 1
     %if binoc has resstarted, reopen pipes
     d = dir('/tmp/binocisnew');
     if isfield(d,'datenum')
         if d.datenum > lastts || isempty(lastts)
             lastts = now;
             go = 1;
         elseif (now-d.datenum) > 1/(24 * 60) % more than a minute with no binco
%shoudl only get here if verg gets no response from binoc. If that is true and 
%/tmp/binocisnew is old, it suggests binoc has died.  But double check
%before autorestarting
             if DATA.autorestart
                 [a, pstr] = system('ps -e | grep binoclean');
                 if isempty(strfind(pstr, 'binoclean.app'))
                     fprintf('Binoc is Not Running - ');
                     if DATA.autoreopen == 0
                         DATA.autoreopen = 2;
                     end
                     myprintf(DATA.frombinocfid,'-show','Restarting Binoc at %s',datestr(now));
                     [a,b] = system('open /local/bin/binoclean.app');
                 end
             end
         end
         
     end
else
    go = 2;
end
if go && DATA.newbinoc ~= 2
     if DATA.autoreopen
         fprintf('Reopening pipes\n');
         if go ==1 %if received NewBinoc, no need to pause
             pause(1);
         end
         DATA = ReadIO(DATA,[],6);
         if isfield(DATA,'reopenstr') && ~isempty(DATA.reopenstr)
             outprintf(DATA,'%s\n',DATA.reopenstr);
         end
         if DATA.autoreopen ==2 %once only
             DATA.autoreopen = 0;
         end
     else
         fprintf('Binoc Restarted - reopen pipes\n');
     end
end
     
function OpenUffFile(a,b, type)
    DATA = GetDataFromFig(a);
    outprintf(DATA,'\\openuff\n');
    myprintf(DATA.cmdfid,'#File %s\n',GetValue(DATA,'uf'));


function Expt = ExptSummary(DATA)
    Expt.first = DATA.Trial.Trial;
    Expt.Stimvals.et = DATA.exptype{1};
    Expt.Stimvals.e2 = DATA.exptype{2};
    Expt.Stimvals.e3 = DATA.exptype{3};
    Expt.Stimvals.ei = DATA.binoc{1}.ei;
    Expt.Stimvals.i2 = DATA.binoc{1}.i2;
    Expt.Stimvals.i3 = DATA.binoc{1}.i3;
    Expt.Stimvals.st = DATA.stimulusnames{DATA.stimtype(1)};
    Expt.Stimvals.Bs = DATA.stimulusnames{DATA.stimtype(2)};
    Expt.Start = now;
    Expt.Header.Name = GetName(DATA.binoc{1}.uf,'-silent');
    Expt.Header.Subject = GetValue(DATA,'monkey');

        
function DATA = RunButton(a,b, type)
        DATA = GetDataFromFig(a);
        if isstruct(a)
            if ~isempty(b)
                caller = b;
            else
                caller = ['Seq'];
            end
        else
            caller = get(a,'String');
            if sum(strcmp(caller,{'Cancel' 'Finish'})) %really hit button
                DATA.canceltimer = 0;
            end
        end
        ts = now;
        myprintf(DATA.frombinocfid,'-show','%s Hit Inexpt %d, type %d %s\n',caller,DATA.inexpt,type,datestr(now));
        myprintf(DATA.cmdfid,'# %s Hit Inexpt %d, type %d %s\n',caller,DATA.inexpt,type,datestr(now));
        PauseRead(DATA,1);
        pauseread = getappdata(DATA.toplevel, 'PauseReading');
        if pauseread ==0
            myprintf(DATA.frombinocfid,'PAUSE ERROR pause %d\n',pauseread);
        end

        DATA.newexptdef = 0;
        DATA.matexpres = [];
        forcestop = sum(strcmp(caller,'ForceCancel'));
        inexpt = 0;
        if type == 1
            if DATA.inexpt == 0 && ~forcestop %sarting a new one. Increment counter
                DATA = CheckPenLog(DATA);

%if storage is off, and correct smr file is open, check that user means it                
                if isfield(DATA.binoc{1},'uf') && ~isempty(DATA.smrfile)
                    if strncmp(DATA.smrfile,DATA.binoc{1}.uf,length(DATA.binoc{1}.uf))
                        if DATA.optionflags.ts == 0
                            yn = gui.Dlg(sprintf('smr file %s is open - Did you want storage on?',DATA.smrfile),DATA.toplevel);
                            if strcmp(yn,'Yes')
                                DATA.optionflags.ts = 1;
                                SendCode(DATA,'optionflag');
                                SetGui(DATA);
                            end
                        end
                    end
                end
%if storage is on, but filname does not match, warn used                
                if DATA.optionflags.ts == 1
                    if ~strncmp(DATA.smrfile,DATA.binoc{1}.uf,length(DATA.binoc{1}.uf))
                        yn = gui.Dlg(sprintf('smr file is %s. Proceed?',DATA.smrfile),DATA.toplevel);
                        if strcmp(yn,'No')
                            return;
                        end                        
                    end
                end
                pauseread = getappdata(DATA.toplevel, 'PauseReading');
                if pauseread ==0
                    myprintf(DATA.frombinocfid,'PAUSE ERROR pause %d\n',pauseread);
                end
                if DATA.optionflags.exm && ~isempty(DATA.matexpt)
                    DATA.matexpres = binoceval(DATA, DATA.matexpt);
                    if isfield(DATA.matexpres,'abort') && DATA.matexpres.abort > 0 %matlab script finds a problem
                        vergwarning(sprintf('%s Says abort',DATA.matexpt));
                        PauseRead(DATA,0);
                        return;
                    end
                    SendManualExpt(DATA);
                end
                DATA = SendManualVals(DATA,'All');
                outprintf(DATA,'\\expt\n');
                DATA.nexpts = DATA.nexpts+1;
                DATA.Expts{DATA.nexpts} = ExptSummary(DATA);
                DATA.optionflags.do = 1;
                DATA.exptstoppedbyuser = 0;
                DATA.trialcounts(6) = 0;
                caller = 'Run';
                inexpt = 1;
                %            DATA = GetState(DATA);
            else
                inexpt = 0;
                DATA.rptexpts = 0;
                if DATA.verbose(4)
                    fprintf('Before Cancel: Inexpt %d\n',DATA.inexpt);
                end
                DATA.exptstoppedbyuser = 2;
                outprintf(DATA,'\\ecancel\n');
            end
        elseif type == 2
            DATA.exptstoppedbyuser = 1;
            outprintf(DATA,'\\estop\n');
            inexpt = 0;
        end
%if expt is over, EXPTOVER should be received. - query state then
%    DATA = GetState(DATA);
%    DATA = ReadFromBinoc(DATA);
    if DATA.verbose(4)
        fprintf('RunButton: Inexpt %d\n',DATA.inexpt);
    end
%    DATA = GetState(DATA,'Runbutton');
    x(1) = mytoc(ts);
    ts = now;
%Waiting for expt end is special case. Binoc may return info before the final
%trial ends, and state is still in expt. So wait for binoc to report that
%expt state has changed.
    while DATA.inexpt ~= inexpt  && mytoc(ts) < 10 %not caught up yet
        pauseread = getappdata(DATA.toplevel, 'PauseReading');
        if pauseread ==0
            myprintf(DATA.frombinocfid,'PAUSE ERROR pause %d\n',pauseread);
        end
        [DATA, str] = ReadFromBinoc(DATA,'RunButton','paused');
        myprintf(DATA.frombinocfid,'Read Done %.2f %s. DATA.inexpt %d\n',mytoc(ts),datestr(now),DATA.inexpt);
        if DATA.inexpt ~= inexpt
            pause(0.5);
        end
    end
    x(2) = mytoc(ts);

    if sum(strcmp(caller,{'Cancel' 'Seqcancel' 'ForceCancel'}))
        if DATA.nexpts <1
            DATA.nexpts = 1;
        end
          DATA = AddTextToGui(DATA,'Cancelled','norec');
                if DATA.nexpts > 0
                DATA.Expts{DATA.nexpts}.last = DATA.Trial.Trial;
                DATA.Expts{DATA.nexpts}.End = now;
                end
                DATA.optionflags.do = 0;
                DATA.exptstoppedbyuser = 1;
                if strcmp(caller,'Cancel')
                    DATA.exptstoppedbyuser = 1;
                end
    elseif strcmp(caller,'Run')
       CheckExptIsGood(DATA);
    elseif strcmp(caller,'Finish')
        if DATA.nexpts <1
            DATA.nexpts = 1;
        end
            DATA.rptexpts = 0;
            DATA.Expts{DATA.nexpts}.last = DATA.Trial.Trial;
            DATA.Expts{DATA.nexpts}.End = now;
            ti = 1 + DATA.Trial.Trial-DATA.Expts{DATA.nexpts}.first;
            DATA.optionflags.do = 0;
            DATA = AddTextToGui(DATA,['Stopped after ' num2str(ti) ' Trials']);
            DATA.exptstoppedbyuser = 1;
    end
    if DATA.inexpt ~= inexpt
            myprintf(DATA.frombinocfid,'-show','Expt Button (%s) Confused Inexpt is %d,%d. Trying again\n',caller,DATA.inexpt,inexpt);
            DATA = GetState(DATA, 'RunButton');
            if DATA.inexpt ~= inexpt
                DATA.confused = 1;
            end
    else
        DATA.confused = 0;
    end
    set(DATA.toplevel,'UserData',DATA);    
    pauseread = getappdata(DATA.toplevel, 'PauseReading');
    if pauseread == 0
        myprintf(DATA.frombinocfid,'PAUSE ERROR pause %d\n',pauseread);
    end
    myprintf(DATA.frombinocfid,'Runbutton Done %.2f pause %d\n',x,pauseread);
    SetGui(DATA);
    PauseRead(DATA,0);

    CheckTimer(DATA);
     
    function TestIO(a,b, str)
        
        onoff = {'off' 'on'};
        DATA = GetDataFromFig(a);
        if ~isfield(DATA,'outpipe')
            DATA = OpenPipes(DATA, 1);
        end
        if strcmp(str,'cancel')
            DATA.canceltimer = 10;
        elseif strcmp(str,'freeze')
            outprintf(DATA,'!freeze');
            if ~isfield(DATA,'freezestimulus')
                DATA.freezestimulus =1;
            else
                DATA.freezestimulus = ~DATA.freezestimulus;
            end
            set(a,'checked',onoff{DATA.freezestimulus+1});
            hm = findobj(allchild(DATA.toplevel),'flat','type','uimenu','Tag','StepFrame');
            if isempty(hm)
               uimenu(DATA.toplevel,'Label','Step','Tag','StepFrame','callback',{@TestIO,'step'});
            end
        elseif strcmp(str,'step')
            fprintf('Stepping frame\n');
            outprintf(DATA,'!step\n');
        elseif strcmp(str,'saveimage')
            outprintf(DATA,'!saveimage\n');
        end
        set(DATA.toplevel,'UserData',DATA);
        
function ElectrodeMoved(F, pos, varargin) %called by ServoDrive
%ServoDrive sends position in microns. Binoc wants mm            
if ~isfigure(F) %if verg was close before servowindow
    return;
end
DATA = get(F, 'UserData');
if strcmp(pos,'close') %Servo Contoller Closing
    DATA.servofig = 0;
    if isfield(varargin{1},'alldepths')
        X = CheckDepthData(varargin{1});
        if isfield(X,'newtimes') && X.newtimes(1) > max(X.alltimes)
            DATA.servodata.alldepths = cat(2,DATA.servodata.alldepths,X.alldepths,X.newdepths);
            DATA.servodata.alltimes = cat(2,DATA.servodata.alltimes,X.alltimes,X.newtimes);
        else
            DATA.servodata.alldepths = cat(2,DATA.servodata.alldepths,X.alldepths);
            DATA.servodata.alltimes = cat(2,DATA.servodata.alltimes,X.alltimes);
        end
        DATA.servodata = CheckDepthData(DATA.servodata);
        DATA.servodata.stepsize = X.stepsize;
    end
    if isfield(DATA,'servotimer')
        stop(DATA.servotimer);
    end
else
    try
    if ~isempty(varargin) && isfield(varargin{1},'alldepths')
        X = CheckDepthData(varargin{1});
        if isempty(DATA.servodata.alltimes)
            DATA.servodata.alltimes = X.alltimes;
            DATA.servodata.alldepths = X.alldepths;
        else
            id = find(X.alltimes > max(DATA.servodata.alltimes));
            DATA.servodata.alldepths = cat(2,DATA.servodata.alldepths,X.alldepths(id));
            DATA.servodata.alltimes = cat(2,DATA.servodata.alltimes,X.alltimes(id));        
        end
    end
    outprintf(DATA,'!seted=%.3f\n',pos./1000);
    DATA.binoc{1}.ed = pos./1000;
    S = DATA.binoc{1};
    if isfield(S,'Pn') && S.Pn >0
        ServoDrive('label', sprintf('Penetration %d at %.1f,%.1f',S.Pn,S.Xp,S.Yp));
    else
        ServoDrive('label', sprintf('Penetration Not Set',S.Pn,S.Xp,S.Yp));
    end
    catch ME
        CheckExceptions(ME);
    end
end
set(DATA.toplevel,'UserData',DATA);



function DATA = CheckDepthData(DATA)
if length(DATA.alltimes) > length(DATA.alldepths)
    DATA.alltimes = DATA.alltimes(1:length(DATA.alldepths));
end


function ElectrodePopup(a,b, fcn, varargin)
  DATA = GetDataFromFig(a);

  if ~isempty(DATA.servoport)
      args = {};
      
%If DATA.servodata.alldepths exists, then Verg has already communicated with 
%ServoDrive.  Might still want to check that it has not been turned off and
%lost position
      if ~isempty(DATA.servodata.alldepths)
          args = {'depths' DATA.servodata.alldepths DATA.servodata.alltimes};
      else
          args = {'startdepth' DATA.binoc{1}.ed .* 1000};
      end
      if ~isempty(DATA.servodata.stepsize)
          args = {args{:} 'stepsize' DATA.servodata.stepsize};
      end
      X = ServoDrive('ttyname',DATA.servoport,'callback',{@ElectrodeMoved, DATA.toplevel},args{:});
      DATA.servofig = X.toplevel;
      figure(DATA.servofig); %force to front
      
      hm = findobj(allchild(DATA.servofig),'flat','type','uimenu','Tag','MarkMenu');
      if isempty(hm)
          hm = uimenu(gcf,'label','Mark','Tag','MarkMenu');
          AddMarkMenu(hm);
      end

      DATA.servotimer = timer('timerfcn',{@CheckServo, DATA.toplevel},'Period',1,'executionmode','fixedspacing');
      start(DATA.servotimer);
      set(DATA.toplevel,'UserData',DATA);
  end
  

function CheckExptIsGood(DATA)
        
function DATA = CheckPenLog(DATA)
%if saving DATA, and servodrive has been moved, Then check penlog is open    
    S = DATA.binoc{1};
    if DATA.optionflags.ts && length(DATA.servodata.alldepths) > 1 %check pen log is open
        if ~isfield(S,'Pn') || S.Pn <= 0
            if DATA.state.penwarned < 2
                vergwarning('Penetration log is not open.');
            else
                cprintf('red','WARNING!!       Penetration log not opened\n');
            end
            DATA.state.penwarned = DATA.state.penwarned+1;
            set(DATA.toplevel,'UserData',DATA);
        else %log is ok. enable warning if this stops being true
            DATA.state.penwarned = 0;
        end
    end
    
    
    
    
function PenLogPopup(a,b)
  DATA = GetDataFromFig(a);
  cntrl_box = findobj('Tag',DATA.windownames{9},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      if strcmp(b,'show')
          if isfield(DATA.binoc{1},'Pn') && DATA.binoc{1}.Pn > 0
              set(cntrl_box,'name',sprintf('Penetration %d',DATA.binoc{1}.Pn));
          else
              set(cntrl_box,'name','Penetration Log (No log open)');
          end
      end
      return;
  end
if length(DATA.winpos{9}) ~= 4
    DATA.winpos{9} = get(DATA.toplevel,'position');
end
cntrl_box = figure('Position', DATA.winpos{9},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{9},'Name','Penetration Log','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
            set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);
    set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);

    nr = 5;
    nc = 7;
    bp = [0.01 0.99-1/nr 1./nc 0.98./nr];

    

    
    bp(1) = 0.01;
    bp(3) = 0.3;
    uicontrol(gcf,'style','text','string','Penetration Number', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.98./nc;
    
    uicontrol(gcf,'style','edit','string',num2str(GetValue(DATA,'Pn')), ...
        'units', 'norm', 'position',bp,'Tag','Pn','callback',{@TextGui, 'Pn'});

    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.1;
    uicontrol(gcf,'style','text','string','X', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.98./nc;
    uicontrol(gcf,'style','edit','string',num2str(GetValue(DATA,'Xp')), ...
        'units', 'norm', 'position',bp,'Tag','Xp','callback',{@TextGui, 'Xp'});
    
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.1;
    uicontrol(gcf,'style','text','string','Y', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.98./nc;
    uicontrol(gcf,'style','edit','string',num2str(GetValue(DATA,'Yp')), ...
        'units', 'norm', 'position',bp,'Tag','Yp','callback',{@TextGui, 'Yp'});

    bp(1) = 0.01;
    bp(2) = bp(2)-1./nr;
    uicontrol(gcf,'style','text','string','Electrode', ...
        'units', 'norm', 'position',bp);

    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 3./nc;
   uicontrol(gcf,'style','pop','string',DATA.electrodestrings, ...
        'units', 'norm', 'position',bp,'value',DATA.electrodeid,'Tag','ElectrodeType','callback',{@MenuGui});
 
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.2;
    uicontrol(gcf,'style','text','string','Impedance', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 1./nc;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.eZ), ...
        'units', 'norm', 'position',bp,'Tag','Impedance','callback',{@MenuGui});
    

    
    bp(2) = bp(2)-1./nr;
    bp(1) = 0.01;
    bp(3) = 0.9./nc;
    uicontrol(gcf,'style','text','string','User', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','userlabel');
    bp(1) = 1./nc;
    bp(3) = 1.9./nc;
    id = find(strcmp(DATA.binoc{1}.ui,DATA.userstrings));
    if isempty(id)
        id = 1;
    end
    uicontrol(gcf,'style','pop','string',DATA.userstrings, ...
        'units', 'norm', 'position',bp,'value',id,'Tag','Experimenter','callback',{@MenuGui});

    bp(1) = 3./nc;
    bp(3) = 0.9./nc;
    uicontrol(gcf,'style','text','string','Monkey', ...
        'units', 'norm', 'position',bp);

    bp(1) = bp(1)+bp(3);
    bp(3) = 1.9./nc;
    id = find(strcmp(DATA.binoc{1}.monkey,DATA.monkeystrs));
    if isempty(id)
        id = 1;
    end
    uicontrol(gcf,'style','pop','string',DATA.monkeystrings, ...
        'units', 'norm', 'position',bp,'value',id,'Tag','Monkey','callback',{@MenuGui});

    bp(2) = bp(2)-1./nr;
    bp(1) = 0.01;
    bp(3) = 0.2;
    uicontrol(gcf,'style','text','string','Tube Protrusion', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.1;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.ePr), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','ePr','callback',{@MenuGui});

    bp(1) = bp(1)+bp(3);
    bp(3) = 0.2;
    uicontrol(gcf,'style','text','string','Coarse mm', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.1;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.coarsemm), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','coarsemm','callback',{@MenuGui});

    bp(1) = bp(1)+bp(3);
    bp(3) = 0.2;
    uicontrol(gcf,'style','text','string','Adapter', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.2;
    uicontrol(gcf,'style','edit','string',DATA.binoc{1}.adapter, ...
        'units', 'norm', 'position',bp,'value',1,'Tag','adapter','callback',{@MenuGui});

    bp(2) = bp(2)-1./nr;
    bp(1) = 0.01;
    bp(3) = 0.15;
    
    
    id = find(strcmp(DATA.binoc{1}.hemi,{'Left' 'Right' 'NotSet'}));
    if isempty(id)
        id = 1;
    end
    
    uicontrol(gcf,'style','text','string','Hemisphere', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.15;
    uicontrol(gcf,'style','pop','string','Left|Right|Unknown', ...
        'units', 'norm', 'position',bp,'value',id,'Tag','hemisphere','callback',{@MenuGui});


    bp(1) = bp(1)+bp(3);
        uicontrol(gcf,'style','text','string','Area', ...
        'units', 'norm', 'position',bp);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.2;
    uicontrol(gcf,'style','pop','string','V1|V2|MT|Vc (calcarine V1)|Unknown', ...
        'units', 'norm', 'position',bp,'value',5,'Tag','VisualArea','callback',{@MenuGui});

    
    bp(1) = 0.99-bp(3);
    uicontrol(gcf,'style','pushbutton','string','Apply', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','Penset','callback',@OpenPenLog);
    bp(1) = bp(1)-bp(3)-0.01;
    uicontrol(gcf,'style','pushbutton','string','Plot', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','PlotPen','callback',@OpenPenLog);
    hm = uimenu(gcf,'label','Mark');
    AddMarkMenu(hm);
    hm = uimenu(gcf,'label','Set');
    uimenu(hm,'label','New Penetration','tag','NewPen','callback',{@MenuBarGui});
    uimenu(hm,'label','Show Penetration Log','tag','ShowPen','callback',{@ShowPenLog, 'popup'});
   
    
set(gcf,'CloseRequestFcn',{@CloseWindow, 9});

function MarkComment(a,b,txt);
  DATA = GetDataFromFig(a);
  if isfield(DATA,'callingfigure')
      DATA = get(DATA.callingfigure,'UserData');
  end
  AddComment(DATA,txt, 'penmenu');
  outprintf(DATA,'cm=%s\n',txt);
        
function AddMarkMenu(hm)
    uimenu(hm,'label','Head Restrained','callback',{@MarkComment 'Head Restrained'});
    uimenu(hm,'label','Entered Brain','callback',{@MarkComment 'Entered Brain'});
    uimenu(hm,'label','Entered GM','callback',{@MarkComment 'GM'});
    uimenu(hm,'label','Entered WM','callback',{@MarkComment 'WM'});
    uimenu(hm,'label','Penetration Missed Lunate/Calcarine','callback',{@MarkComment 'MissedDeepSulci'});
        
function OptionPopup(a,b)
  DATA = GetDataFromFig(a);
  cntrl_box = findobj('Tag',DATA.windownames{2},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      return;
  end
if length(DATA.winpos{2}) ~= 4
    DATA.winpos{2} = get(DATA.toplevel,'position');
end
f = fields(DATA.optionflags); %Don't mess with order of thes
silentf = fields(DATA.silentoption);
nc = 4;
nr = ceil((length(f)+nc-length(silentf))/nc);
scrsz = get(0,'Screensize');
cntrl_box = figure('Position', DATA.winpos{2},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{2},'Name','Options','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
        set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
            set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);

bp = [0.01 0.99-1/nr 1./nc 1./nr];
nf = 0;
for j = 1:length(f)
    if ~isfield(DATA.silentoption,f{j})
        nf = nf+1;
        bp(1) = floor((nf-1)/nr) .* 1./nc;
        bp(2) = 1- ((1+rem(nf-1,nr)) .* 1./nr);
        if strncmp('lbl',f{j},3)
            uicontrol(gcf,'style','text','string',[DATA.optionstrings.(f{j}) ':'], ...
                'units', 'norm', 'position',bp,'value',DATA.optionflags.(f{j}),'Tag',f{j},...
                'callback',{@HitToggle, f{j}},'foregroundcolor','r');
        else
            uicontrol(gcf,'style','checkbox','string',[DATA.optionstrings.(f{j}) '(' f{j} ')'], ...
                'units', 'norm', 'position',bp,'value',DATA.optionflags.(f{j}),'Tag',f{j},'callback',{@HitToggle, f{j}});
        end
    end
end
f = fields(DATA.stimflags{1});
for j = 1:length(f)
    str = DATA.stimflagnames.(f{j});
    k = nf+j; 
    bp(1) = floor(k/nr) .* 1./nc;
    bp(2) = 1- (rem(k,nr) .* 1./nr);
    uicontrol(gcf,'style','checkbox','string',str, ...
        'units', 'norm', 'position',bp,'value',DATA.stimflags{1}.(f{j}),'Tag',f{j},'callback',{@StimToggle, f{j}});

end

function CommentPopup(a,b,type)
  DATA = GetDataFromFig(a);
  src = a;
  
  wn = find(strncmp(DATA.windownames,'comment',6));  
  cntrl_box = findobj('Tag',DATA.windownames{wn},'type','figure');
  if ~isempty(cntrl_box);
      lst = findobj(cntrl_box,'tag','CommentList');
  end
  if strncmp(type,'show',4)
      f = strrep(type,'show','');
      DATA.show.comment.(f) = ~DATA.show.comment.(f);
      if DATA.show.comment.(f)
          set(a,'checked','on');
      else
          set(a,'checked','off');
      end
      SetData(DATA);
      CommentPopup(lst, [], 'update')
  elseif strcmp(type,'update')
      if ~isempty(cntrl_box);
          strs = {};
          for j = 1:length(DATA.Comments)
              strs{j} = DATA.Comments(j).comment;
              if DATA.show.comment.time
                  strs{j} = [datestr(DATA.Comments(j).date,'hh:mm') ': ' strs{j}];
              end                      
              if DATA.show.comment.depth
                  strs{j} = sprintf('ed%.3f: %s', DATA.Comments(j).ed, strs{j});
              end
          end
          if ~isempty(strs)
              set(lst,'string',strs);
          end
      end
  end
  if ~strncmp(type,'popup',5)
      return;
  end
  
  cntrl_box = findobj('Tag',DATA.windownames{wn},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      return;
  end
  cntrl_box = figure('Position', DATA.winpos{wn},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{wn},'Name','Comments','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
    set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
    set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);

    hm = uimenu(cntrl_box,'Label','Show');
    sm = uimenu(hm,'Label','Comments','callback',{@CommentPopup, 'showcomments'},'checked','on');
    sm = uimenu(hm,'Label','Time','callback',{@CommentPopup, 'showtime'},'checked','on');
    sm = uimenu(hm,'Label','Electrode Depth','callback',{@CommentPopup, 'showdepth'},'checked','on');

    lst = uicontrol(gcf, 'Style','list','String', 'Code List :* = more help with mouse click',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
        'Tag','CommentList',...
        'callback',{@CommentPopup, 'help'},...
        'units','norm', 'Position',[0.01 0.085 0.99 0.91]);
    
if ~isfield(DATA,'show') || ~isfield(DATA.show,'comment')    
    DATA.show.comment.time = 1;
    DATA.show.comment.depth = 1;
    DATA.show.comment.comment = 1;
    SetData(DATA);
end
CommentPopup(lst, [], 'update')

function ShowPenLog(a,b,type)
  DATA = GetDataFromFig(a);
  src = a;
  
  wn = find(strncmp(DATA.windownames,'showpen',6));  
  cntrl_box = findobj('Tag',DATA.windownames{wn},'type','figure');
  if ~isempty(cntrl_box);
      lst = findobj(cntrl_box,'tag','PenetrationLog');
  end
  if strncmp(type,'show',4)
      f = strrep(type,'show','');
      DATA.show.penlog.(f) = ~DATA.show.penlog.(f);
      if DATA.show.penlog.(f)
          set(a,'checked','on');
      else
          set(a,'checked','off');
      end
      SetData(DATA);
      ShowPenLog(lst, [], 'update')
  elseif strcmp(type,'update')
      if ~isempty(cntrl_box);
          txt = scanlines(['/local/' DATA.binoc{1}.monkey '/pen' num2str(DATA.binoc{1}.Pn) '.log']);
          
          strs = {};
          for j = 1:length(txt)
              if length(txt{j}) > 0
                  go = 1;
              else
                  go = 0;
              end
              if strncmp(txt{j},'ed',2) && DATA.show.penlog.depth == 0
                  go = 0;
              elseif strncmp(txt{j},'Rewards',6) && DATA.show.penlog.rewards == 0
                  go = 0;
              elseif strncmp(txt{j},'Expt',4) && DATA.show.penlog.expts == 0
                  go = 0;
              elseif strfind(txt{j},'bwticks') 
                  go = 0;
              elseif ~isempty(strfind(txt{j},' File ')) && DATA.show.penlog.otherlines == 0
                  go = 0;
              elseif strncmp(txt{j},'cm=rf',5) && DATA.show.penlog.autocomments == 0
                  go = 0;
              elseif strncmp(txt{j},'cm=rf',5) && DATA.show.penlog.autocomments == 0
                  go = 0;
              elseif sum(strncmp(txt{j},{'Experimenter' 'Hemisphere' 'VisualArea' 'Electrode'},8)) ...
                 && DATA.show.penlog.otherlines == 0
             go = 0;
              end
              if go
                  strs{end+1} = txt{j};
              end
          end
          if ~isempty(strs)
              set(lst,'string',strs);
          end
      end
  end
  if ~strncmp(type,'popup',5)
      return;
  end
  
  cntrl_box = findobj('Tag',DATA.windownames{wn},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      return;
  end
  cntrl_box = figure('Position', DATA.winpos{wn},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{wn},'Name','Penetration Log Contents','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
    set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
    set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);

    if ~isfield(DATA,'show') || ~isfield(DATA.show,'penlog')
        DATA.show.penlog.depth= 1;
        DATA.show.penlog.rewards = 1;
        DATA.show.penlog.expts= 1;
        DATA.show.penlog.autocomments= 0;
        DATA.show.penlog.otherlines= 1;
        SetData(DATA);
    end
    onoff = {'off' 'on'};
    hm = uimenu(cntrl_box,'Label','Show');
    sm = uimenu(hm,'Label','Electrode Depth','callback',{@ShowPenLog, 'showdepth'},'checked',onoff{DATA.show.penlog.depth+1});
    sm = uimenu(hm,'Label','Rewards','callback',{@ShowPenLog, 'showrewards'},'checked',onoff{DATA.show.penlog.rewards+1});
    sm = uimenu(hm,'Label','Expts','callback',{@ShowPenLog, 'showexpts'},'checked',onoff{DATA.show.penlog.expts+1});
    sm = uimenu(hm,'Label','Auto Comments','callback',{@ShowPenLog, 'showautocomments'},'checked',onoff{DATA.show.penlog.autocomments+1});
    sm = uimenu(hm,'Label','Other Lines','callback',{@ShowPenLog, 'showotherlines'},'checked',onoff{DATA.show.penlog.otherlines+1});

    lst = uicontrol(gcf, 'Style','list','String', 'Penetration Log',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
        'Tag','PenetrationLog',...
        'callback',{@ShowPenLog, 'help'},...
        'units','norm', 'Position',[0.01 0.085 0.99 0.91]);
    

ShowPenLog(lst, [], 'update')


function CodesPopup(a,b, type)  

  DATA = GetDataFromFig(a);
  src = a;
  e = {};
  if isnumeric(type) | strmatch(type,{'bycode' 'bylabel' 'bygroup' 'numeric' 'printcodes' 'byhelp'},'exact')
      lst = findobj(get(get(a,'parent'),'parent'),'Tag','CodeListString');
      if isnumeric(type)
          if type == 8
            set(lst,'string','Codes controlling PsychoPhysics/Reward/Eye movement');
          else
          end
          id = find(bitand([DATA.comcodes.group],type) > 0);
        [c,b] = sort({DATA.comcodes(id).code});
        b = id(b);
      elseif strcmp(type,'printcodes')
          F = get(a,'parent');
         [outname, path] = uiputfile([DATA.localmatdir '/BinocCodes.txt'], 'Save Binoc Codes');
         if ~ischar(outname) %user cancelled
             return;
         end
         outname = [path outname];
         it = findobj(F,'tag','CodeListString');
         txt = get(it,'String');
         fid = fopen(outname,'w');
         for j = 1:size(txt,1)
             fprintf(fid,'%s\n',deblank(txt(j,:)));
             code = regexprep(txt(j,:),'\s.*','');
             if isfield(DATA.helpkeys.extras,code)
                 for k = 1:length(DATA.helpkeys.extras.(code))
                     fprintf(fid,'     %s\n',DATA.helpkeys.extras.(code){k});
                 end
             end
         end
         fclose(fid);
         fprintf('Codes writted to %s\n',outname);
         return;
      elseif strcmp(type,'bycode')
          set(lst,'string','Alphabetical by code.  :* = more help with mouse click');
          [c,b] = sort({DATA.comcodes.code});
          e = setdiff(fields(DATA.helpstrs),{DATA.comcodes.code}); %help on things not in comcodes
      elseif strcmp(type,'bylabel')
          set(lst,'string','Alphabetical by Label :* = more help with mouse click');
          [c,b] = sort({DATA.comcodes.label});
          e = setdiff(fields(DATA.helpstrs),{DATA.comcodes.code}); %help on things not in comcodes
      elseif strcmp(type,'byhelp')
          set(lst,'string','Alphabetical by Help :* = more help with mouse click');
          f = fields(DATA.helpstrs);
          for j = 1:length(f)
              helpstr{j} = DATA.helpstrs.(f{j});
              hid = find(strcmp(f{j},{DATA.comcodes.code}));
              if ~isempty(hid)
                  cid(j) = hid;
              else
                  cid(j) = NaN;
              end
          end
          [c, b] = sort(helpstr);
          helpcodes = f(b);
          b = cid(b);
      elseif strcmp(type,'bygroup')
          set(lst ,'string','Groups: :* = more help with mouse click');
          id = find(bitand(1,[DATA.comcodes.group]));
          [c,b] = sort({DATA.comcodes(id).code});
          labels{1} = 'Stimulus Rendering';
          id = find(bitand(2,[DATA.comcodes.group]));
          [c,d] = sort({DATA.comcodes(id).code});
          labels{2} = 'Psychophysics/Reward';
          b = [b 0 d];
          id = find(bitand(4,[DATA.comcodes.group]));
          [c,d] = sort({DATA.comcodes(id).code});
          labels{3} = 'Group 3';
          b = [b 0 d];
          id = find(bitand(8,[DATA.comcodes.group]));
          [c,d] = sort({DATA.comcodes(id).code});
          labels{4} = 'Group 4';
          b = [b 0 d];
          id = find(bitand(16,[DATA.comcodes.group]));
          [c,d] = sort({DATA.comcodes(id).code});
          labels{5} = 'Group 5';
          b = [b 0  d];
      else
          set(lst,'string','Numerical');
          b = 1:length(DATA.comcodes);
      end
      nl=1;
      a = get(lst,'string');
      nlab = 0;
      nc = 0;
      for j = 1:length(b)
          if b(j) > 0
              code = DATA.comcodes(b(j)).code;
              codetype = DATA.comcodes(b(j)).group;
          elseif isnan(b(j))
              code = helpcodes{j};
              codetype = 0;
          else
              code = '';
              codetype = 0;
          end
          if ~strcmp(code,'xx')
              ns = max([5 - length(code) 1]);
              ns = 1+ round(ns-1) .* 1.6;
              nc = nc+1;
              if b(j) > 0
                  s = sprintf('%s%*s%s',code,ns,' ',DATA.comcodes(b(j)).label);
              elseif isnan(b(j))  %help with no code
                  f = helpcodes{j};
                  s = sprintf('%s %s',code,DATA.helpstrs.(f));
              else
                  nlab = nlab+1;
                  a(nc+nl,1) = ' ';
                  nl = nl+1;
                  s = labels{nlab};
              end
              if isfield(DATA.helpstrs,code)
                  if isfield(DATA.helpkeys.extras,code)
                      s = regexprep(s,' ',' *','once');
                  end
                  s = [s '   ;   ' DATA.helpstrs.(code)];
                  keys{j+nl} = DATA.helpkeys.KeyWords.(code);
              end
              a(nc+nl,1:length(s)) = s;
              a(nc+nl,2+length(s):end) = 0;
          end
      end
      for j = 1:length(e)
          nc = nc+1;
          code = e{j};
          if DATA.helpstrs.(code)(1) == '!'
            s = ['!' code '   ;   ' DATA.helpstrs.(code)(2:end)];
          else
            s = [code '   ;   ' DATA.helpstrs.(code)];
          end
          a(nc+nl,1:length(s)) = s;
          a(nc+nl,2+length(s):end) = 0;
      end
      a = a(1:nc+nl,:);
      cmenu = uicontextmenu;
      set(cmenu,'UserData',lst);
      uimenu(cmenu,'label','Help','Callback',{@CodeListMenu, 'Help'});
      a = AddOptionHelp(DATA,a);
      set(lst,'string',a);
      setappdata(GetFigure(lst),'HelpKeyList',keys);
      setappdata(GetFigure(lst),'OldText',a);
      %set(lst,'uicontextmenu',cmenu);    
      
  end
  if ~strncmp(type,'popup',5)
      return;
  end
  cntrl_box = findobj('Tag',DATA.windownames{4},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      return;
  end
  if DATA.togglecodesreceived == 0 %no binoc connection
      X = load('comcodes.mat');
      DATA.comcodes = X.comcodes;
  end
  
  cntrl_box = figure('Position', DATA.winpos{4},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{4},'Name','Code list','menubar','none');
    
    if isfield(DATA,'toplevel')
        setappdata(cntrl_box,'ParentFigure',DATA.toplevel);
    else
          DATA.toplevel = cntrl_box;
          SetData(DATA);
    
    end
        set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
        set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);

    hm = uimenu(cntrl_box,'Label','List by');
    sm = uimenu(hm,'Label','By Code','callback',{@CodesPopup, 'bycode'});
    sm = uimenu(hm,'Label','By Label','callback',{@CodesPopup, 'bylabel'});
    sm = uimenu(hm,'Label','By Group','callback',{@CodesPopup, 'bygroup'});
    sm = uimenu(hm,'Label','Psych/Reward','callback',{@CodesPopup, 8 });
    sm = uimenu(hm,'Label','Numerical','callback',{@CodesPopup, 'numeric'});
    sm = uimenu(hm,'Label','By Help','callback',{@CodesPopup, 'byhelp'});
    helpmenu = sm;
    hm = uimenu(cntrl_box,'Label','Print','callback',{@CodesPopup, 'printcodes'});
    hm = uimenu(cntrl_box,'Label','Search');
    sm = uimenu(hm,'Label','Ignore case','callback',{@SearchList, 'IgnoreCase'},'Tag','IgnoreCase','checked','on');
    
    uicontrol(gcf,'style','pop','string','Search: All|Search: Codes|Search: Labels|Search: Help|Search: HelpFile','tag','SearchMode',...
        'units','norm', 'Position',[0 0.01 0.3 0.08]);

    srch = uicontrol(gcf, 'Style','edit','String', '',...
        'callback', @SearchList, ...
        'units','norm', 'Position',[0.3 0.01 0.99 0.08]);

    lst = uicontrol(gcf, 'Style','list','String', 'Code List :* = more help with mouse click',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
        'Tag','CodeListString',...
        'callback',@HelpHit,...
        'units','norm', 'Position',[0.01 0.085 0.99 0.91]);
a = get(lst,'string');
nc = 1;
keys = {};
for j = 1:length(DATA.comcodes)
    code = DATA.comcodes(j).code;
    if ~strcmp(code,'xx')
        ns = max([5 - length(code) 1]);
        ns = 1+ round(ns-1) .* 1.6;
        s = sprintf('%s%*s%s',code,ns,' ',DATA.comcodes(j).label);
        if isfield(DATA.helpstrs,code)
            s = [s '  : ' DATA.helpstrs.(code)];
            keys{j} = DATA.helpkeys.KeyWords.(code);
        end
        nc = nc+1;
        a(nc,1:length(s)) = s;
    else
        nc = nc;
    end
end
a = AddOptionHelp(DATA,a);
a = AddVergHelp(DATA,a);

set(lst,'string',a);
setappdata(GetFigure(lst),'HelpKeyList',keys);
setappdata(GetFigure(lst),'OldText',a);

if strcmp(type,'popuphelp')
    CodesPopup(helpmenu,b,'byhelp');
elseif strcmp(type,'popup')
    CodesPopup(helpmenu,b,'bycode');
end

function HelpHit(a,b)
DATA = GetDataFromFig(a);

str = get(a,'string');
l = get(a,'value');
if length(l) == 1
    code = str(l,:);
    code = regexprep(code,'\s.*','');
    fpos = get(GetFigure(a),'position');
    if isfield(DATA.helpkeys.extras,code)
        cm = uicontextmenu;
        X = DATA.helpkeys.extras.(code);
        for j = 1:length(X)
            uimenu(cm,'label',X{j});
        end
        set(cm,'position', [50 fpos(4)/2],'visible','on');
    end
end


function a = AddVergHelp(DATA,a)
    
   f = setdiff(fields(DATA.helpstrs),{DATA.comcodes.code});
   for j = 1:length(f)
       s = sprintf('%s %s\n',f{j},DATA.helpstrs.(f{j}));
       a(end+1,1:length(s)) = s;
   end

function a = AddOptionHelp(DATA,a)
        
if isfield(DATA.helpkeys,'options')
    s = 'Binary options (op=+xx to set. op=-xx to unset)';
    a(end+1,1:length(s)) = s;
    
    f = fields(DATA.optionflags);
    for j = 1:length(f)
        if isfield(DATA.helpkeys.options,f{j})
            s = sprintf('+%s %s   ;   %s',f{j},DATA.optionstrings.(f{j}),DATA.helpkeys.options.(f{j}));
        else
            s = sprintf('+%s %s\n',f{j},DATA.optionstrings.(f{j}));
        end
        a(end+1,1:length(s)) = s;
    end

         
    f = fields(DATA.helpkeys.options);
    for j = 1:length(f)
        s = sprintf('+%s %s',f{j},DATA.helpkeys.options.(f{j}));
        a(end+1,1:length(s)) = s;
    end
end

function ToggleCheck(a)
    if strcmp(get(a,'checked'),'off')
        set(a,'checked','on');
    else
        set(a,'checked','off');
    end

function SearchList(a,b, varargin)

    
    tag = get(a,'tag');
    if strcmp(tag,'IgnoreCase')
        ToggleCheck(a);
        return;
    end
    newpattn = 0;
    pttn = get(a,'string');
    lastpttn = get(a,'UserData');
    set(a,'UserData',pttn);    
    F = get(a,'parent');
    findn = getappdata(F,'findn');
    it = findobj(F,'tag','CodeListString');
    if isempty(findn) %hit return again
        findn = 1;
    elseif ~strcmp(pttn,lastpttn)%new serach
        newpattn = 1;
    else
        setappdata(F,'findn',[]);
        if isappdata(F,'OldText')
            txt = getappdata(F,'OldText');
            set(it,'string',txt);
            set(F,'Name','Code List');
            return;
        end
    end
    if strcmp(pttn, lastpttn)
        findn = findn+1;
    end
    
    ignorecase = 0;
    mit = findobj(F,'tag','IgnoreCase');
    if strcmp(get(mit,'checked'),'on')
        ignorecase = 1;
    end
    sm = findobj(F,'tag','SearchMode');
    searchmode = get(sm,'value');
    if newpattn
        txt = getappdata(F,'OldText');
    else
        txt = get(it,'String');
        setappdata(F,'OldText',txt);
    end
    keys = getappdata(F,'HelpKeyList');
    keys{size(txt,1)+1} = '';
    found = [];
    for j = 1:size(txt,1)
        t = txt(j,:);
        if isempty(strfind(t,':'))
            t = [t ' : '];
        end
        if searchmode == 2 %just codes
            s = regexprep(t,'\s.*','');
        elseif searchmode == 3 %just labels
            s = regexprep(t,'\s(.*):.*','$1');
        elseif searchmode == 4 %just help
            s = regexprep(t,'.* : ','');
        else
            s = [txt(j,:) keys{j}];
        end
        if ignorecase
            x = strfind(lower(s),lower(pttn));
        else
            x = strfind(s,pttn);
        end
        if ~isempty(x)
            found(j) = 1;
        end
    end
    found = find(found);
    if ~isempty(found)
       if findn > length(found)
           findn = 1;
       end
       if length(found) == 0
%           set(it,'value',found(findn));
       else
           str = sprintf('%d matches for %s',length(found),pttn);
           txt(1,1:length(str)) = str;
           txt(1,(length(str)+1):end) = ' ';
           set(it,'value',1);
           set(it,'string',txt([1 found],:));
       end
        set(F,'Name',sprintf('%d Matches for %s',length(found),pttn));
    else
        set(F,'Name',sprintf('No Matches for %s',pttn));
    end
    setappdata(F,'findn',findn);
    
function StatusPopup(a,b, type)  

  DATA = GetDataFromFig(a);
  if ~strcmp(type,'popup')
      return;
  end
  cntrl_box = findobj('Tag',DATA.windownames{5},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      return;
  end
  cntrl_box = figure('Position', DATA.winpos{5},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{5},'Name','Status Lines from binoc','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
        set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);

    lst = uicontrol(cntrl_box, 'Style','list','String', 'Code List',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
         'Tag','NextButton',...
'units','norm', 'Position',[0.01 0.01 0.99 0.99]);
set(lst,'string',DATA.Statuslines,'fontsize',DATA.font.FontSize, 'FontName', DATA.font.FontName);
DATA.statusitem = lst;
if ~isfield(DATA,'showstatus') || ~isfield(DATA.showstatus,'trials')
    DATA.showstatus.trials = 1;
    DATA.showstatus.status = 1;
    DATA.showstatus.errors = 1;    
    DATA.showstatus.expt = 1;
    DATA.showstatus.comment = 1;
    DATA.showstatus.update = 1;
end
f = fields(DATA.showstatus)
onoff = {'off' 'on'};
cm = uimenu(cntrl_box,'Label','Show');
for j = 1:length(f)
    uimenu(cm,'label',f{j},'tag',f{j},'checked',onoff{1+DATA.showstatus.(f{j})},'callback',@SetShowStatus);
end
set(DATA.toplevel,'UserData',DATA);

function SetShowStatus(a,b)
DATA = GetDataFromFig(a);
tag = get(a,'tag');
onoff = {'off' 'on'};
DATA.showstatus.(tag) = ~DATA.showstatus.(tag);
set(a,'checked',onoff{1+DATA.showstatus.(tag)});
ShowStatusStrings(DATA);
SetData(DATA);

function ShowStatusStrings(DATA)
    if ishandle(DATA.statusitem)
        f = {'status' 'trials' 'errors' 'expt' 'comment'};
        for j = 1:length(f)
            showlines(j) = DATA.showstatus.(f{j});
        end
        showlines = find(showlines);
        id = find(ismember(DATA.statustypes,showlines));
        set(DATA.statusitem,'string',DATA.Statuslines(id),'value',length(id));
    end
        
function DATA = RunExptSequence(DATA, str, line)

    nread = 0;
    DATA.matlabwasrun = 0;
    if ~iscellstr(str)
        str = cellstr(str);
    end
    runlines = find(strncmp('!expt',str,5) | strncmp('!trial',str,5)) ;
    if isempty(runlines)
        warndlg('Sequence must contain !expt or !trial','Sequence file error');
        lastline = 1;
    else
        lastline = runlines(end);
    end
    firstline = line;
    if DATA.runsequence == 0
        return;
    end
 
    if line > lastline
        if DATA.rptexpts > 0
            DATA.rptexpts = DATA.rptexpts-1;
            fprintf('Repeating sequence, %d to go\n',DATA.rptexpts);
            if DATA.restartbinoc 
                PauseRead(DATA,1);
                DATA= RestartBinoc(DATA);
                %force a short pause to communicate with new binoc=
                fprintf('Pausing for Restart of sequence %s\n',datestr(now));
                DATA = uipause(now, max([1 DATA.binoc{1}.seqpause]),'Fixed Sequence Pause', DATA);
                PauseRead(DATA,0);
                firstline = 1; %don't have second pause
            end            
            line = 1;
        else
            DATA.seqline = 0;
            myprintf(DATA.cmdfid,'-show','Sequence End\n');
            return;
        end
    end
    if ischar(str)
        astr = str;
        clear str;
        for j = 1:size(astr,1)
            str{j} = deblank(astr(j,:));
        end
    end
    lastline = '';
for j = line:length(str)
    if DATA.verbose(4)
        fprintf('From Seq window %s\n',str{j});
    end
    nread = 1+j-line;
    rptstogo = DATA.rptexpts;
    DATA = InterpretLine(DATA,str{j},'fromseq');
    DATA = uipause(DATA.pausetime,DATA.readpause,'Pause in sequence', DATA); %for pauses set in window
    DATA.readpause = 0;
    DATA.rptexpts = rptstogo; %don't let stimfiles change this
    if strncmp(str{j},'!expt',5)
        if DATA.inexpt
            cprintf('red','!!Sequence !expt at line %d but already in Expt!!!',j);
            return;
        end
%need to do this before sending !expt to binoc, so that UserData is set
% before binoc calls back with settings
        cntrl_box = findobj('Tag',DATA.windownames{8},'type','figure');
          lst = findobj(cntrl_box,'Tag','SequenceList');
%        set(lst,'value',j);
        if strcmp(str{j},'!exptrpt') %repeat exact expt, without calling matexpt again
            DATA.matlabwasrun = 1;
        end
        if firstline > 1
            fprintf('Pausing for Next Expt in Sequnce Line %d %s\n',j,datestr(now));
            DATA= uipause(now, DATA.binoc{1}.seqpause,'Fixed Sequence Pause', DATA);
        end
  %if mat was called in the sequence file, don't want it overridden by the matept file
        if DATA.optionflags.exm && ~isempty(DATA.matexpt) && DATA.matlabwasrun == 0
            fprintf('Running %s\n',DATA.matexpt);
            DATA.matexpres = [];
            DATA.matexpres = binoceval(DATA, DATA.matexpt);            
            DATA.matlabwasrun = 1;
            if isfield(DATA.matexpres,'abort') && DATA.matexpres.abort > 0
                vergwarning(sprintf('%s Says abort',DATA.matexpt));
                return;
            end
            SendManualExpt(DATA);
        end
        myprintf(DATA.cmdfid,'!expt line %d',j);
        DATA.nexpts = DATA.nexpts+1;
        DATA.Expts{DATA.nexpts} = ExptSummary(DATA);
        DATA.seqline = j;
        if DATA.rptexpts > 0
            xs = sprintf(' %d rpts to go',DATA.rptexpts);
        else
            xs = '';
        end
            
        AddTextToGui(DATA,sprintf('#Expt %s Line %d of Sequence%s',Expt2Name(DATA.Expts{DATA.nexpts}),j,xs));

        if j > 1 && j <= length(str)
            s = str{j-1};
        elseif j > 0
            s = str{1};
        else
            s = 'end';
        end
        set(cntrl_box,'Name',sprintf('Sequence at line %d (%s)',DATA.seqline,s));
        set(DATA.toplevel,'UserData',DATA);
        outprintf(DATA,'#From RunSequence\n');
        outprintf(DATA,'%s\n',str{j}); %this runs expt in binoc
        DATA = AddStatusLine(DATA,sprintf('RunSequence Line %d. at %s %d Repeats left',j,datestr(now,'hh:mm.ss'),DATA.rptexpts),4);
        return;
    elseif strncmp(str{j},'!trial',5)
        outprintf(DATA,'!onetrial\n');
        pause(0.2);
        DATA = DrainBinocPipe(DATA,'waitforstim');
    elseif strncmp(str{j},'!mat',4)
        DATA = AddStatusLine(DATA, sprintf('RunSequence Line %d: %s',str{j}),4);
        if isfield(DATA.matexpres,'abort') && DATA.matexpres.abort > 0
            vergwarning(sprintf('%s Says abort',DATA.matexpt));
            return;
        end
    end
    if ~sum(strncmp(str{j},'expt',4)) %don't send these lines to binoc
        outprintf(DATA,'%s#RunSeq\n',str{j});
    end
    DATA = LogCommand(DATA, str{j},'norec');
    lastline = str{j};
end

function SendManualExpt(DATA)
    
    SendCode(DATA,'exp');
    X = DATA.matexpres;
    if isfield(X,'binocstrs')
        for j = 1:length(X.binocstrs)
            outprintf(DATA,'%s\n',X.binocstrs{j});
        end
    end
        
    if isfield(X,'exvals')
        a = unique(X.exvals(:,1));
        DATA.binoc{1}.nt = length(a);
        SendCode(DATA,'nt');
        for j = 1:length(a)
            outprintf(DATA,'EA%d=%s\n',j-1,num2str(a(j)));
        end
        if size(X.exvals,2) > 1
            a = unique(X.exvals(:,2));
            DATA.nstim(2) = length(a);
            SendCode(DATA,'n2');
            for j = 1:length(a)
                outprintf(DATA,'EB%d=%s\n',j-1,num2str(a(j)));
            end
        end
        if size(X.exvals,2) > 2
            a = unique(X.exvals(:,3));
            DATA.nstim(3) = length(a);
            SendCode(DATA,'n3');
            for j = 1:length(a)
                outprintf(DATA,'EC%d=%s\n',j-1,num2str(a(j)));
            end
        end
    end
    outprintf(DATA,'EDONE\n');

function DATA = uipause(start, secs, msg, DATA)

    if secs <= 0
        return;
    end
    args = {'frompause'};
    paused = getappdata(DATA.toplevel,'PauseReading');
    if paused
        args = {args{:} 'paused'};
    end
    days = secs/(24 * 60 * 60);
    wh = waitbar(0,sprintf('%.1f sec %s at %s',secs,msg,datestr(start)));
    while (now - start) < days
        dt = (now - start)/days;
        wdur = (now-start) * (24 * 60 * 60);
        waitbar(dt,wh,sprintf('%.1f/%.1f sec %s at %s',wdur,secs,msg,datestr(start)));
        if nargin > 3
            DATA = ReadFromBinoc(DATA, args{:});
        end
    end
    delete(wh);


function DATA = ContinueSequence(DATA, varargin)
    showlog = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'log',3)
            showlog = 1;
        end
        j = j+1;
    end
  cntrl_box = findobj('Tag',DATA.windownames{8},'type','figure');
  lst = findobj(cntrl_box,'Tag','SequenceList');
  if DATA.newbinoc == 0 && ~isempty(lst) %don't call this when just parsing initial state
      if showlog
          str = get(lst,'string');
          if ~iscell(str)   %!!! Seems not to be consistent
              str = cellstr(str);
          end
          if length(str) <= DATA.seqline
              str{DATA.seqline+1} = 'Back to start';
          end
          myprintf(DATA.cmdfid,'#Sequence continuing from line %d\n#%s\n',DATA.seqline,str{DATA.seqline+1});
      end
      DATA = RunExptSequence(DATA,get(lst,'string'),DATA.seqline+1);
  end

function DATA = LogCommand(DATA, str, varargin)
    j = 1;
    reccmd = 1;
    d = [' #' datestr(now,'HH:MM')];
    while  j <= length(varargin)
        if strcmp(varargin{j},'norec')
            reccmd = 0;
        elseif strcmp(varargin{j},'notime')
            d = '';
        end
        j = j+1;
    end

    if reccmd
        DATA.commands = {DATA.commands{:} str};
        DATA.commandctr = length(DATA.commands)+1;
    end
    
    if DATA.cmdfid > 0
        fprintf(DATA.cmdfid,'%s%s\n',str,d);
    end
    if nargout == 0
        set(DATA.toplevel,'UserData',DATA);
    end
    
    
function name = Seq2Name(str)
    name = '??';
    id = find(strncmp('!mat',str,4));
    if ~isempty(id)
        x = str{id(1)};
        name = regexprep(x(6:end),'(.*','');
    end
    
    
function SequencePopup(a,exptlines,type, varargin)

  DATA = GetDataFromFig(a);
  cntrl_box = findobj('Tag',DATA.windownames{8},'type','figure');
  if ~strncmp(type,'popup',5)
      if strcmp(type,'run')
          DATA = GetState(DATA,'runseq');
          DATA.runsequence = 1;
          DATA.seqline = 0;
          lst = findobj(cntrl_box,'Tag','SequenceList');
          DATA = RunExptSequence(DATA,get(lst,'string'),1);
      elseif strcmp(type,'Revert')
          oldid = varargin{1};
          lst = findobj(cntrl_box,'Tag','SequenceList');
          set(lst,'string',DATA.sequences(oldid).str);
      elseif strcmp(type,'Save')
          lst = findobj(cntrl_box,'Tag','SequenceList');
          str = get(lst,'string');
          if ~isfield(DATA,'sequences')
              it= uimenu(gcf,'Label','Revert','Tag','SequenceRevert');
              n = 1;
          else
              n = length(DATA.sequences)+1;
              it = findobj(allchild(cntrl_box),'flat','Tag','SequenceRevert');
          end
          DATA.sequences(n).str = str;
          DATA.sequences(n).name = Seq2Name(str);
          DATA.sequences(n).time = now;
          if ~isempty(it)
              uimenu(it,'Label',sprintf('%s %s',DATA.sequences(n).name,datestr(now,'HH:MM')),'callback',{@SequencePopup,'Revert', n});
          end
          SetData(DATA);  
      elseif strcmp(type,'pause')
          str = get(a,'String');
          if strcmp(str,'Pause');
              set(a,'String','continue');
              DATA.runsequence = 0;
          else
              set(a,'String','Pause');
              DATA.runsequence = 1;
%if continue is hit while still running the same expt where pause was hit,
%dont advance to the next experiment immediatel - wait for expt to finish
              if DATA.inexpt == 0
                  DATA = ContinueSequence(DATA,'log');
              end
          end
          SetData(DATA);
      end
      return;
  end
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      if strcmp(type,'popupnew')
          lst = findobj(allchild(cntrl_box),'flat','tag', 'SequenceList');
          SetSequenceString(lst, DATA, exptlines);
      end
      return;
  end
  
  nr=10; nc=4;
  cntrl_box = figure('Position', DATA.winpos{8},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{8},'Name','Expt Sequence','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
        set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);

        bp(1) = 0.01;
bp(2) = 1-2./nr;
bp(3) = 1./nc;
bp(4) = 2./nr;
uicontrol(gcf,'style','pushbutton','string','Run', ...
    'Callback', {@SequencePopup, 'run'} ,...
    'units', 'norm', 'position',bp,'value',1);

if Dev(DATA)
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','Pause', ...
        'Callback', {@SequencePopup, 'pause'} ,...
        'units', 'norm', 'position',bp,'value',1);
   
end
bp(1) = bp(1)+bp(3)+0.01;
uicontrol(gcf,'style','pushbutton','string','Save', ...
    'Callback', {@SequencePopup, 'Save'} ,...
    'units', 'norm', 'position',bp,'value',1);

if isfield(DATA,'sequences')
    it= uimenu(gcf,'Label','Revert','Tag','SequenceRevert');
    for j = 1:length(DATA.sequences)
        uimenu(it,'Label',sprintf('%s %s',DATA.sequences(j).name,datestr(DATA.sequences(j).time,'HH:MM')),'callback',{@SequencePopup,'Revert', j});
    end
end

usejava =0;
if usejava
lst = jcontrol(gcf,'javax.swing.JTextField',...
                    'Units','normalized',...
             'Tag','SequenceList',...
                    'Position',[0.01 0.01 0.98 0.99-2./nr]);
      set(lst,'Text','this is a test');
else
    lst = uicontrol(gcf, 'Style','edit','String', 'sequence',...
        'HorizontalAlignment','left',...
        'Max',10,'Min',0,...
        'Tag','SequenceList',...
        'units','norm', 'Position',[0.01 0.01 0.99 0.99-2./nr]);
    SetSequenceString(lst,DATA, exptlines);
end
set(DATA.toplevel,'UserData',DATA);

function SetSequenceString(lst, DATA, exptlines)
    if iscellstr(exptlines) && ~isempty(exptlines)
        set(lst,'string',exptlines);
    elseif strncmp(class(exptlines),'matlab.ui',9) || isempty(exptlines)
        sid = find(strcmp('sequence',DATA.exptlines));
        if ~isempty(sid)
            set(lst,'string',DATA.exptlines(sid+1:end));
        elseif isfield(DATA,'matexpt')
            exptline{1} = sprintf('!mat=%s',DATA.matexpt);
            exptline{2} = '!expt';
            set(lst,'string',exptline);
        end
    end

function d = Dev(DATA)
    if isfield(DATA.binoc{1},'ui') && strcmp(DATA.binoc{1}.ui,'bgc')
        d = true;
    else
        d = false;
    end
        
function CodeListMenu(a,b,c)
    lst = get(get(a,'parent'),'userdata');
    strs = get(lst,'string')
    line = get(lst,'value');
    code = regexprep(strs(line,:),' .*','');
    DATA = GetDataFromFig(lst);
    if ~isfield(DATA,'helpdata')
        DATA.helpdata = loadhelp('/local/binochelp.txt');
    end
    id = find(strcmp(code,DATA.helpdata.codes));
    if ~isempty(id)
       helpdlg(DATA.helpdata.help{id},code); 
    else
       helpdlg(sprintf('No Help on %s :(',code),code); 
    end
    
function helpdata = loadhelp(name)
    
    fid = fopen(name,'r');
    if (fid > 0)
        a = textscan(fid,'%s','delimiter','\n');
    end
    txt = a{1};
    hlp = [];
    nh = 0;
    for j = 1:length(txt)
        if txt{j}(1) == ':'
            id = strfind(txt{j},':');
            if nh > 0
                helpdata.codes{nh} = code;
                helpdata.help{nh} = hlp;
            end
            nh = nh+1;
            code = txt{j}(id(1)+1:id(2)-1);
            hlp = txt{j}(id(2):end);
            
        else 
            hlp = [hlp txt{j}];
        end
        helpdata.codes{nh} = code;
        helpdata.help{nh} = hlp;
    end
        
    
 function DATA = UpdateLogFile(DATA)
  
     if ~isfield(DATA.binoc{1},'lo') || isempty(DATA.binoc{1}.lo)
         return;
     end
     name = DATA.binoc{1}.lo;
     if isfield(DATA,'Trials')
         if length(DATA.Trials) < 10
             return;
         end
     end
    fid = fopen(name,'a');
    if fid < 0
        fprintf('Cannot Append to %s\n',name);
        return;
    end
    first = DATA.Trials(1).Start;
    last = DATA.Trials(end).Start;
    dur = (last-first)./24;
    fprintf(fid,'%s %d Trials, %.2f hrs %s to %s Rw %.1f. User %s/%s\n',datestr(now),length(DATA.Trials),dur,datestr(first),datestr(last),DATA.binoc{1}.Trw,DATA.binoc{1}.ui,GetUserName());
    fclose(fid);


function DATA = ReadLogFile(DATA, name)

    fid = fopen(name,'r');
    if fid < 0
        fprintf('Cannot Read %s\n',name);
        return;
    end
    s = textscan(fid,'%s','delimiter','\n');
    s = s{1};
    gotwt = 0;
    savetime = 0;
    for j = 1:length(s)
        e = strfind(s{j},'=');
        if ~isempty(e)
            value = s{j}(e+1:end);
            code = s{j}(1:e-1);
        else 
            code = [];
        end
        if strncmp(s{j},'Saved',5)
            %used to read back to last save, but that assumes it was all
            %there
            savetime = datenum(s{j}(8:end));
        elseif strncmp(s{j},'Gain',4)
            DATA.Coil.gain = sscanf(s{j}(6:end),'%f');
        elseif strncmp(s{j},'Offset',6)
            DATA.Coil.offset = sscanf(s{j}(8:end),'%f');
        elseif strncmp(s{j},'Phase',5)
            DATA.Coil.phase = sscanf(s{j}(7:end),'%f');
        elseif strncmp(s{j},'CriticalWeight',14)
            DATA.Coil.CriticalWeight = sscanf(s{j}(15:end),'%f');
            DATA.Coil.CriticalDate = savetime;
        elseif strncmp(s{j},'so',2)
            DATA.Coil.so = sscanf(s{j}(4:end),'%f');
        elseif strncmp(s{j},'we',2)
            we = sscanf(s{j}(3:end),'%f');
            if we(1) > 0
                gotwt = gotwt+1;
                DATA.binoc{1}.we = we(1);
                DATA.Coil.weightdate = savetime;
                DATA.Coil.lastwt = we(1);
            end
            if length(we) > 1 && we(2) > 0
                DATA.Coil.CriticalWeight = we(2);
                DATA.Coil.CriticalDate = savetime;
            end
        elseif sum(strncmp(s{j},{'MicroDrive'},6)) 
            if isempty(code) %no =, use first space as delimiter
                code = regexprep(s{j},'\s.*','');
                value = strrep(s{j},code,'');
                code = genvarname(code);
            end
            if ~isempty(code)
                DATA.Coil.Xtra.(code) = value;
            end
        end
    end
    fclose(fid);
    if gotwt
        SendCode(DATA,'we');
    end
    
function MonkeyLogPopup(a,b, type, channel)
  DATA = GetDataFromFig(a);
  

  if DATA.Coil.gain(1) == 0 %not read yet
      DATA = ReadLogFile(DATA, DATA.binoc{1}.lo);
  end
  
  if ~strcmp(type,'popup')
  value = str2num(get(a,'string'));

  id = strmatch(type,{'RH' 'LH' 'RV' 'LV' 'save' 'clear' 'popup'});
  if strncmp(type,'Offset',6)
      DATA.Coil.offset(channel) = value;
  elseif strncmp(type,'getsoft',5)
      DATA.Coil.so = DATA.binoc{1}.so;
     SetTextItem(gcf,'RH',DATA.binoc{1}.so(1));
     SetTextItem(gcf,'LH',DATA.binoc{1}.so(2));
     SetTextItem(gcf,'RV',DATA.binoc{1}.so(3));
     SetTextItem(gcf,'LV',DATA.binoc{1}.so(4));
  elseif strncmp(type,'applysoft',8)
      DATA.binoc{1}.so = DATA.Coil.so;
      SendCode(DATA,'so');
  elseif strncmp(type,'CriticalWeight',8)
      DATA.Coil.CriticalWeight = value;
      DATA.Coil.CriticalDate = now;
  elseif strncmp(type,'Comment',8)
      if ishandle(a)
          DATA.Coil.Comment = get(a,'string');
      end
  elseif strncmp(type,'Weight',5)
      DATA.binoc{1}.we = value;
      DATA.Coil.Weight = value;
      DATA.Coil.weightdate = now;
  elseif strncmp(type,'phase',5)
      DATA.Coil.phase(channel) = value;
  elseif strncmp(type,'Gain',4)
      DATA.Coil.gain(channel) = value;
  elseif strncmp(type,'offset',6)
      DATA.Coil.offset(channel) = value;
  elseif strncmp(type,'soft',4)
      DATA.Coil.so(channel) = value;
  elseif strcmp(type,'savelog')
      fid = fopen(DATA.binoc{1}.lo,'a');
      if fid > 0
      fprintf(fid,'Saved: %s\n',datestr(now));
      fprintf(fid,'so%s\n',sprintf(' %.2f',DATA.Coil.so));
      fprintf(fid,'Gain%s\n',sprintf(' %.2f',DATA.Coil.gain));
      fprintf(fid,'Offset%s\n',sprintf(' %.2f',DATA.Coil.offset));
      fprintf(fid,'Phase%s\n',sprintf(' %.2f',DATA.Coil.phase));
      if  strcmp(type,'savelog') && DATA.Coil.Weight > 0 % only save weight if it has been added in GUI.
          fprintf(fid,'we%.2f',DATA.Coil.Weight);
          if DATA.Coil.CriticalDate - now > 0.5
              fpritnf(fid,' %.2f\n',DATA.Coil.CriticalWeight);
              fprintf('Setting Critical WEight to %.2f\n',DATA.Coil.CriticalWeight);
          else
              fprintf(fid,'\n');
          end
      else
          fprintf('No Weight Given\n');
      end
      if isfield(DATA.Coil,'Comment') && ~isempty(DATA.Coil.Comment)
          fprintf(fid,'#%s\n',DATA.Coil.Comment);
      end
      fclose(fid);
      fprintf('Saved to %s\n',DATA.binoc{1}.lo);
      end
      close(GetFigure(a));
  end
  set(DATA.toplevel,'UserData',DATA);
      return;
  end
  cntrl_box = findobj('Tag',DATA.windownames{6},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      return;
  end
if length(DATA.winpos{6}) ~= 4
    DATA.winpos{6} = get(DATA.toplevel,'position');
end
cntrl_box = figure('Position', DATA.winpos{6},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{6},'Name','monkeylog','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
    set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
    set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);
    
    if isfield(DATA.Coil,'Xtra')        
        nr = 8;
    else
        nr = 7;
    end
nc=6;

bp = [0.01 0.99-1/nr .99/nc 1./nr];


    bp = [2/nc 0.99-1/nr 0.99./nc 1./nr];
uicontrol(gcf,'style','text','string','RH', ...
        'units', 'norm', 'position',bp,'value',1);
             bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','text','string','LH', ...
        'units', 'norm', 'position',bp,'value',1);


             bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','text','string','RV', ...
        'units', 'norm', 'position',bp,'value',1);

         bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','text','string','LV', ...
        'units', 'norm', 'position',bp,'value',1);

    bp(1) = bp(1)+bp(3)+0.01;
        bp = [0.01 0.99-2/nr 1.99/nc 1./nr];
     uicontrol(gcf,'style','text','string','SoftOff', ...
        'units', 'norm', 'position',bp,'value',1);   

bp = [2/nc 0.99-2/nr 0.99./nc 1./nr];
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.so(1)), ...
        'Callback', {@MonkeyLogPopup, 'soft', 1},'Tag','RH',...
        'units', 'norm', 'position',bp);


bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.so(2)), ...
        'Callback', {@MonkeyLogPopup, 'soft' 2},'Tag','LH',...
        'units', 'norm', 'position',bp);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.so(3)), ...
        'Callback', {@MonkeyLogPopup, 'soft', 3},'Tag','RV',...
        'units', 'norm', 'position',bp);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.so(4)), ...
        'Callback', {@MonkeyLogPopup, 'soft' 4},'Tag','LV',...
        'units', 'norm', 'position',bp);    

        bp = [0.01 0.99-3/nr 1.99/nc 1./nr];
     uicontrol(gcf,'style','text','string','Gain', ...
        'units', 'norm', 'position',bp,'value',1);   

bp = [2/nc 0.99-3/nr 0.99/nc 1./nr];
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.gain(1)), ...
        'Callback', {@MonkeyLogPopup, 'GainRH', 1},'Tag','GainRH',...
        'units', 'norm', 'position',bp);
   
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.gain(2)), ...
        'Callback', {@MonkeyLogPopup, 'GainLH', 2},'Tag','GainLH',...
        'units', 'norm', 'position',bp);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.gain(3)), ...
        'Callback', {@MonkeyLogPopup, 'GainRV', 3},'Tag','GainRV',...
        'units', 'norm', 'position',bp);
    
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.gain(4)), ...
        'Callback', {@MonkeyLogPopup, 'GainLV' 4},'Tag','GainLV',...
        'units', 'norm', 'position',bp);    

    bp = [0.01 0.99-4/nr 1.99/nc 1./nr];
     uicontrol(gcf,'style','text','string','Phase', ...
        'units', 'norm', 'position',bp,'value',1);   
bp = [2/nc 0.99-4/nr 0.99/nc 1./nr];
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.phase(1)), ...
        'Callback', {@MonkeyLogPopup, 'phaseRH', 1},'Tag','phaseRH',...
        'units', 'norm', 'position',bp);
   
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.phase(2)), ...
        'Callback', {@MonkeyLogPopup, 'phaseLH', 2},'Tag','phaseLH',...
        'units', 'norm', 'position',bp);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.phase(3)), ...
        'Callback', {@MonkeyLogPopup, 'phaseRV', 3},'Tag','phaseRV',...
        'units', 'norm', 'position',bp);
    
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.phase(4)), ...
        'Callback', {@MonkeyLogPopup, 'phaseLV', 4},'Tag','phaseLV',...
        'units', 'norm', 'position',bp);    
    
    
    bp = [0.01 0.99-5/nr 1.99/nc 1./nr];
     uicontrol(gcf,'style','text','string','Offset', ...
        'units', 'norm', 'position',bp,'value',1);   
    bp = [2/nc 0.99-5/nr 0.99/nc 1./nr];
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.offset(1)), ...
        'Callback', {@MonkeyLogPopup, 'offsetRH', 1},'Tag','offsetRH',...
        'units', 'norm', 'position',bp);
    
   
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.offset(2)), ...
        'Callback', {@MonkeyLogPopup, 'offsetLH', 2},'Tag','offsetLH',...
        'units', 'norm', 'position',bp);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.offset(3)), ...
        'Callback', {@MonkeyLogPopup, 'offsetRV', 3},'Tag','offsetRV',...
        'units', 'norm', 'position',bp);
    
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.offset(4)), ...
        'Callback', {@MonkeyLogPopup, 'offsetLV', 4},'Tag','offsetLV',...
        'units', 'norm', 'position',bp);    

bp(2) = bp(2)- 1./nr;
bp(1) = 0.01;

    if GetValue(DATA,'weightdate') >0
        lstr = sprintf('Weight(%s)',datestr(DATA.Coil.weightdate,'dd/mm/yy'));
    else
        lstr = 'Weight';
    end

    it = uicontrol(gcf,'style','text','string',lstr, ...
    'units', 'norm', 'position',bp,'value',1);

%Dont put current weight variarble in here - user should add new one
bp(1) = bp(1)+bp(3)+0.01;
if DATA.Coil.Weight > 0
    str = sprintf('%.2f',DATA.Coil.Weight);    
elseif isfield(DATA.binoc{1},'we') && DATA.binoc{1}.we > 0
    str = sprintf('(%.1f)',DATA.binoc{1}.we);
elseif DATA.Coil.lastwt > 0
    str = sprintf('(%.2f)',DATA.Coil.lastwt);        
else
    str = '0.00';
end


    uicontrol(gcf,'style','edit','string',str, ...
        'Callback', {@MonkeyLogPopup, 'Weight'},'Tag','Weight',...
        'units', 'norm', 'position',bp);

bp(1) = bp(1)+bp(3)+0.01;
    if isfield(DATA.Coil,'CriticalDate') && DATA.Coil.CriticalDate > 0
        lstr = sprintf('CriticalW(%s)',datestr(DATA.Coil.CriticalDate,'mm/yyyy'));
    else
        lstr = 'Critical Weight';
    end
    it = uicontrol(gcf,'style','text','string',lstr, ...
    'units', 'norm', 'position',bp,'value',1);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.Coil.CriticalWeight), ...
        'Callback', {@MonkeyLogPopup, 'CriticalWeight'},'Tag','CriticalWeight',...
        'units', 'norm', 'position',bp);
bp(1) = bp(1)+bp(3)+0.01;
bp(3) = 0.99-bp(1);
    uicontrol(gcf,'style','edit','string','', ...
        'Callback', {@MonkeyLogPopup, 'Comment'},'Tag','Comment',...
        'units', 'norm', 'position',bp);

bp(3) = 1./nc;
    bp(1) = 0.01;
bp(2) = bp(2)- 1./nr;
uicontrol(gcf,'style','pushbutton','string','Save', ...
    'Callback', {@MonkeyLogPopup, 'savelog'} ,...
    'units', 'norm', 'position',bp,'value',1);
bp(1) = bp(1)+bp(3)+0.01;
bp(3) = 2./nc;
    uicontrol(gcf,'style','pushbutton','string','Apply Softoff', ...
        'Callback', {@MonkeyLogPopup, 'applysoft'} ,...
        'units', 'norm', 'position',bp,'value',1);
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','Use Current Softoff', ...
        'Callback', {@MonkeyLogPopup, 'getsoft'} ,...
        'units', 'norm', 'position',bp,'value',1);
if isfield(DATA.Coil,'Xtra')
    bp(2) = bp(2)- 1./nr;
    f = fields(DATA.Coil.Xtra);
    s = '';
    for j = 1:length(f)
        s = sprintf('%s %s=%s',s,f{j},DATA.Coil.Xtra.(f{j}));
    end
    bp(1) = 0.01;
    bp(3) = 0.99;
    it = uicontrol(gcf,'style','text','string',s, ...
    'units', 'norm', 'position',bp,'value',1);bp(1) = bp(1)+bp(3)+0.01;

    
end
    
set(gcf,'CloseRequestFcn',{@CloseWindow, 6});

function SetSoftOffWindow(DATA, a)
    if ~isfigure(a)
        return;
    end
     SetTextItem(a,'RH',DATA.binoc{1}.so(1));
     SetTextItem(a,'LH',DATA.binoc{1}.so(2));
     SetTextItem(a,'RV',DATA.binoc{1}.so(3));
     SetTextItem(a,'LV',DATA.binoc{1}.so(4));

function SoftoffPopup(a,b, type)
  DATA = GetDataFromFig(a);
  
  
  id = strmatch(type,{'RH' 'LH' 'RV' 'LV' 'null' 'clear' 'popup' 'set'});
  if isempty(id)'
      return;
  elseif id < 5
      DATA.binoc{1}.so(id) = str2num(get(a,'string'));
      SendCode(DATA, 'so');
      set(DATA.toplevel,'UserData',DATA);
  elseif id == 5 || id ==8
      if id == 5
     outprintf(DATA,'sonull\n');
     DATA = ReadFromBinoc(DATA);
      end
      SetSoftOffWindow(DATA, gcf);
  elseif id == 6
     outprintf(DATA,'so=0 0 0 0\n');
     SetTextItem(gcf,'RH',0);
     SetTextItem(gcf,'LH',0);
     SetTextItem(gcf,'RV',0);
     SetTextItem(gcf,'LV',0);
  end
  
  if ~strcmp(type,'popup')
      return;
  end
  cntrl_box = findobj('Tag',DATA.windownames{3},'type','figure');
  if ~isempty(cntrl_box)
      figure(cntrl_box);
      return;
  end
if length(DATA.winpos{3}) ~= 4
    DATA.winpos{3} = get(DATA.toplevel,'position');
end
cntrl_box = figure('Position', DATA.winpos{3},...
        'NumberTitle', 'off', 'Tag',DATA.windownames{3},'Name','Softoff','menubar','none');
    set(cntrl_box,'UserData',DATA.toplevel);
        set(cntrl_box,'DefaultUIControlFontSize',DATA.font.FontSize);
        set(cntrl_box,'DefaultUIControlFontName',DATA.font.FontName);

nr = 2;
nc=6
bp = [0.01 0.99-1/nr 0.115 1./nr];
    uicontrol(gcf,'style','text','string','RH', ...
        'units', 'norm', 'position',bp,'value',1);

bp(1) = bp(1)+bp(3)+0.0
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.so(1)), ...
        'Callback', {@SoftoffPopup, 'RH'},'Tag','RH',...
        'units', 'norm', 'position',bp);
   
         bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','text','string','LH', ...
        'units', 'norm', 'position',bp,'value',1);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.so(2)), ...
        'Callback', {@SoftoffPopup, 'LH'},'Tag','LH',...
        'units', 'norm', 'position',bp);
    
         bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','text','string','RV', ...
        'units', 'norm', 'position',bp,'value',1);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.so(3)), ...
        'Callback', {@SoftoffPopup, 'RV'},'Tag','RV',...
        'units', 'norm', 'position',bp);
    
         bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','text','string','LV', ...
        'units', 'norm', 'position',bp,'value',1);

bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.binoc{1}.so(4)), ...
        'Callback', {@SoftoffPopup, 'LV'},'Tag','LV',...
        'units', 'norm', 'position',bp);    

bp(1) = 0.01;
bp(2) = bp(2)- 1./nr;

    uicontrol(gcf,'style','pushbutton','string','Null', ...
        'Callback', {@SoftoffPopup, 'null'} ,...
        'units', 'norm', 'position',bp,'value',1);
    
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','Clear', ...
        'Callback', {@SoftoffPopup, 'clear'} ,...
        'units', 'norm', 'position',bp,'value',1);

bp(3) = 2./nc;
bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','Save to log', ...
        'Callback', {@MonkeyLogPopup, 'save'} ,...
        'units', 'norm', 'position',bp,'value',1);
    
set(gcf,'CloseRequestFcn',{@CloseWindow, 3});


function CloseWindow(a,b,wid)
  DATA = GetDataFromFig(a);
  x = get(a, 'position');
  DATA.winpos{wid} = x;
  if isfigure(DATA.toplevel)
  set(DATA.toplevel,'UserData',DATA);
  end
  delete(a);
        
    
function StepperPopup(a,b)
  DATA = GetDataFromFig(a);
  cntrl_box = findobj('Tag',DATA.tag.stepper,'type','figure');
if ~isempty(cntrl_box)
    figure(cntrl_box);
    return;
end
scrsz = get(0,'Screensize');
cntrl_box = figure('Position', [10 scrsz(4)-480 300 450],...
        'NumberTitle', 'off', 'Tag',DATA.tag.stepper,'Name','Stepper','menubar','none');
    nr = 6;
    bp = [0.01 0.99-1/nr 0.1 1./nr];
    uicontrol(gcf,'style','pushbutton','string','+', ...
        'Callback', {@Stepper, 1, 1},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','-', ...
        'Callback', {@Stepper, -1, 1},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.3;
    uicontrol(gcf,'style','pop','string','10|20|50|100|200', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','StepSize1');
    bp(1) = 0.01;
    bp(2) = bp(2) - bp(4)-0.01;
    uicontrol(gcf,'style','pushbutton','string','+', ...
        'Callback', {@Stepper, 1, 2},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','-', ...
        'Callback', {@Stepper, -1, 2},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    bp(3) = 0.3;
    uicontrol(gcf,'style','edit','string','10', ...
        'units', 'norm', 'position',bp,'value',1,'Tag','StepSize2');
    bp(1) = 0.01;
    bp(2) = bp(2) - bp(4)-0.01;
    uicontrol(gcf,'style','edit','string',num2str(DATA.stepperpos), ...
        'units', 'norm', 'position',bp,'value',1,'Tag','StepperPosition');
    set(cntrl_box,'UserData',DATA.toplevel);
    
function Stepper(a,b, step, type)
    DATA = GetDataFromFig(a);
    sfig = get(a,'parent');
    it = findobj(sfig,'Tag','StepSize1');
    DATA.stepsize(1) = Menu2Val(it);
    it = findobj(sfig,'Tag','StepSize2');
    DATA.stepsize(2) = Text2Val(it);
    if step > 0
        s = sprintf('ed+%.3f\n',DATA.stepsize(type));
    else
        s = sprintf('ed-%.3f\n',DATA.stepsize(type));
    end
    outprintf(DATA,'%s\n',s);
    if DATA.penid > 0
        fprintf(DATA.penid,'%s\n',s);
    end
        
    ReadFromBinoc(DATA);
    
function val = Menu2Val(it)
val = NaN;
if isempty(it)
    return;
end
j = get(it(1),'value');
s = get(it(1),'string');
val = str2num(s(j,:));

function val = Menu2Str(it)
val = '';
if isempty(it)
    return;
end
j = get(it(1),'value');
s = get(it(1),'string');
if iscellstr(s)
    val = s{j};
else
val = s(j,:);
end

function OpenPenLog(a,b, varargin)
    DATA = GetDataFromFig(a);
    F = get(a,'parent');
    btn = get(a,'tag');
    
    
    if strcmp(btn,'Penset')
        DATA.binoc{1}.Xp = Text2Val(findobj(F,'Tag','Xp'));
        DATA.binoc{1}.Yp = Text2Val(findobj(F,'Tag','Yp'));
        DATA.binoc{1}.Pn = Text2Val(findobj(F,'Tag','Pn'));
        DATA.binoc{1}.ePr = Text2Val(findobj(F,'Tag','ePr'));
        DATA.binoc{1}.eZ = Text2Val(findobj(F,'Tag','Impedance'));
        DATA.binoc{1}.adapter = get(findobj(F,'Tag','adapter'),'string');
        DATA.binoc{1}.hemi = Menu2Str(findobj(F,'Tag','hemisphere'));
        DATA.binoc{1}.Vn = Menu2Str(findobj(F,'Tag','VisualArea'));
        DATA.binoc{1}.ui = Menu2Str(findobj(F,'Tag','Experimenter'));
        DATA.binoc{1}.coarsemm = Menu2Str(findobj(F,'Tag','coarsemm'));
        SendCode(DATA,{'Pn' 'Xp' 'Yp' 'ui' 'Electrode' 'adapter' 'eZ' 'ePr' 'hemi' 'coarsemm' 'Vn'});
        outprintf(DATA,'!openpen');
    elseif strcmp(btn,'PlotPen')
        name = sprintf('/local/%s/pen%d.log',DATA.binoc{1}.monkey,DATA.binoc{1}.Pn);
        GetFigure('Pen');
        hold off;
        PlotOnePen(name,'allcomments');
    end
    PenLogPopup(DATA,'show');    
    S = DATA.binoc{1};
    if isfield(S,'Pn') && S.Pn > 0
        str = sprintf('Penetration %d at %.1f,%.1f',S.Pn,S.Xp,S.Yp);
        set(F,'Name',str);
    end
    set(DATA.toplevel,'UserData',DATA);
    
    
    
function val = Text2Val(it)
val = NaN;
if isempty(it)
    return;
end
s = get(it,'string');
val = str2num(s);

        
function GoToggle(a,b)       
    DATA = GetDataFromFig(a);
    go = get(a,'value');
    if go
      outprintf(DATA,'\\go\n');
    else
      outprintf(DATA,'\\stop\n');
    end
    DATA.optionflags.do  = go;
    set(DATA.toplevel,'UserData',DATA);
    
function HitToggle(a,b, flag)       
    DATA = GetDataFromFig(a);
%    flag = get(a,'Tag');
    DATA.optionflags.(flag) = get(a,'value');
%Don't allow Rmonoc and Lmonoc to be checked at the same time    
    if strcmp(flag,'lm') 
        if DATA.optionflags.lm && DATA.optionflags.rm
            DATA.optionflags.rm = 0;
        end
    elseif strcmp(flag,'rm') 
        if DATA.optionflags.lm && DATA.optionflags.rm
            DATA.optionflags.lm = 0;
        end
    elseif strcmp(flag,'ts') %storage on/off confirm before off mid expt 
        if DATA.optionflags.ts == 0 && DATA.inexpt
            yn = confirm('Sure you want to turn storage off mid expt?');            
            if (yn == 0)
                DATA.optionflags.ts = 1;
            end
        end
    end
    s = 'op=';
    f = fields(DATA.optionflags);
    for j = 1:length(f)
        if DATA.optionflags.(f{j})
            s = [s '+' f{j}];
        else
%            s = [s '-' f{j}];
        end
    end
    if DATA.verbose(4)
        fprintf('op=%s\n',s);
    end
    outprintf(DATA,'op=0\n%s\n',s);
    ReadFromBinoc(DATA);
    CheckTimer(DATA);
    SetGui(DATA,'set');
 
   
    
function SendCode(DATA, code)
    
    if iscellstr(code)
        for j = 1:length(code)
            SendCode(DATA,code{j});
        end
            return;
    end
    if strncmp(code,'Electrode',8)
        if length(DATA.binoc{1}.Electrode)
            outprintf(DATA,'Electrode=%s\n',DATA.binoc{1}.Electrode);
        elseif DATA.electrodeid > 0
            outprintf(DATA,'Electrode=%s\n',DATA.electrodestrings{DATA.electrodeid});
        end
    elseif strcmp(code,'ed')
        x = GetValue(DATA,'ed');
        if ~isnan(x)
            outprintf(DATA,'!seted=%.3f\n',x);
        end
    elseif strcmp(code,'exp')
        if isfield(DATA.matexpres,'stimdir')
            outprintf(DATA,'exp=%s\n',DATA.matexpres.stimdir);
        end
    else
        [s, ~, valid] = CodeText(DATA, code);
        if length(s) && valid >= 0 && s(end) ~= '='
            outprintf(DATA,'%s\n',s);
        end
    end

function val = GetOption(DATA,code)    
    if isfield(DATA,'optionflags') && isfield(DATA.optionflags,code)
        val = DATA.optionflags.(code);
    else
        val = 0;
    end

function val = GetValue(DATA,code)    
%Gets value of a code

if isfield(DATA.binoc{1},code)
    val = DATA.binoc{1}.(code);
elseif strcmp(code,'weightdate')
    if ~isfield(DATA.Coil,code) 
        val = NaN;
    else
        val = DATA.Coil.weightdate;
    end
else
    id = find(strcmp(code,DATA.vergonlycodes));
    if isempty(id)
        val = NaN;
    else
        val = DATA.(DATA.vergonlycodes{id});
    end
end
    
function [s, lbl, type, cid] = CodeText(DATA,code, varargin)
s = [];
lbl = [];
type = 0;
id = 0;

cstim = DATA.currentstim;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'back',4)
        cstim = 2;
    elseif strncmpi(varargin{j},'ChoiceU',7)
        cstim = 3;
    elseif strncmpi(varargin{j},'ChoiceD',7)
        cstim = 4;
    end
    j = j+1;
end


cid = find(strcmp(code,{DATA.comcodes.code}));

if strcmp(code,'optionflag')
        s = 'op=';
        f = fields(DATA.optionflags);
        for j = 1:length(f)
            if DATA.optionflags.(f{j})
                s = [s '+' f{j}];
            else
    %            s = [s '-' f{j}];
            end
        end
        s = sprintf('op=0\n%s\n',s);
        s= [s sprintf('%s',StimToggleString(DATA,DATA.currentstim))];
    elseif strcmp(code,'nr')
        s = sprintf('%s=%d',code,DATA.binoc{1}.nr);
    elseif strmatch(code,{'nt' 'n2' 'n3'})
        id = strmatch(code,{'nt' 'n2' 'n3'});
        s = sprintf('%s=%d',code,DATA.nstim(id));
    elseif strmatch(code,{'et' 'e2' 'e3'})
        id = strmatch(code,{'et' 'e2' 'e3'});
        s = sprintf('%s=%s',code,DATA.exptype{id});
    elseif strcmp(code,'expts')
        s = sprintf('et=%s\nei=%s\nem=%.6f\nnt=%d\n',DATA.exptype{1},DATA.binoc{1}.ei,DATA.mean(1),DATA.binoc{1}.nt);
        s = [s sprintf('e2=%s\ni2=%s\nm2=%6f\nn2=%d\n',DATA.exptype{2},DATA.binoc{1}.i2,DATA.mean(2),DATA.nstim(2))];
        s = [s sprintf('e3=%s\ni3=%s\nm3=%.6f\nn3=%d\n',DATA.exptype{3},DATA.binoc{1}.i3,DATA.mean(3),DATA.nstim(3))];
        s = AddCustomStim(DATA,s,[1:3]);
    elseif strcmp(code,'verbose')
        s = sprintf('verbose=%d\n',DATA.verbose(2));
    elseif strcmp(code,'st')
        s = sprintf('st=%s',DATA.stimulusnames{DATA.stimtype(cstim)});
    elseif strcmp(code,'vve')
        s = sprintf('vve=%s',DATA.vergversion);
    elseif strcmp(code,'pf')
        s = 'pf=';
        f = fields(DATA.showflags);
        s = sprintf('pf=%s',sprintf('+%s',f{:}));
    elseif isfield(DATA.binoc{cstim},code)
        if length(cid) ==1
            lbl = DATA.comcodes(cid).label;
            type = DATA.comcodes(cid).group;
        end
        if length(cid) ==1 && DATA.comcodes(cid).type == 'C'
            s = sprintf('%s=%s',code,DATA.binoc{cstim}.(code));
            if isempty(DATA.binoc{cstim}.(code))
                type = -1;
            end
        else
            s = sprintf('%s=%s',code,num2str(DATA.binoc{cstim}.(code)'));
        end  
    elseif sum(strcmp(code,{DATA.strcodes.code})) == 1
        id = find(strcmp(code,{DATA.strcodes.code}));
        s = sprintf('%s=%s',code,DATA.binocstr.(code));
        lbl = DATA.strcodes(id).label;
    end
        
function StimToggle(a,b, flag)       
    DATA = GetDataFromFig(a);
%    flag = get(a,'Tag');
    DATA.stimflags{1}.(flag) = get(a,'value');
    outprintf(DATA,'%s\n',StimToggleString(DATA,DATA.currentstim));
    ReadFromBinoc(DATA);

function s = StimToggleString(DATA, current)
    s = 'fl=';
    if current > length(DATA.stimflags)
        return;
    end
    f = fields(DATA.stimflags{current});
%always send + and -, so that don't have to track clearing
    for j = 1:length(f)
        if DATA.stimflags{current}.(f{j})
            s = [s '+' f{j}];
        else
            s = [s '-' f{j}];
        end
    end
        
function OtherToggles(a,b,flag)

    DATA = GetDataFromFig(a);
    v= get(a,'value');
    if v
        c = '+';
    else
      c = '-';
    end
    if strcmp(flag,'XYR')
        DATA.showxy(1) = v;
        outprintf(DATA,'ch10%c\n',c);
    elseif strcmp(flag,'XYL')
        DATA.showxy(2) = v;
        outprintf(DATA,'ch11%c\n',c);
    elseif strcmp(flag,'XYB')
        DATA.showxy(3) = v;
        outprintf(DATA,'ch12%c\n',c);
    end        
    set(DATA.toplevel,'UserData',DATA);
    
    
function [DATA, txt] = PrevCommand(DATA, src, step)
    txt = '';
    
    if DATA.newchar ==1  %if typed something, use this to complete commands with
        DATA.completestr = src.Text;
    end
    
    s = get(DATA.txtrec,'string');
    nlines = size(s,1);
    if DATA.commandctr >= length(DATA.commands) && step == 1 && DATA.historyctr > length(DATA.oldcmds)
        DATA.commandctr = length(DATA.commands)+1;
        txt = '';
    elseif (DATA.commandctr > 1 && step == -1) || (DATA.commandctr < length(DATA.commands) && step == 1 && DATA.historyctr == length(DATA.oldcmds)+1)
        if isempty(DATA.completestr)
            DATA.commandctr = DATA.commandctr+step;
        else
            %                fprintf('History with %s\n',DATA.completestr);
            id = find(strncmp(DATA.completestr,DATA.commands,length(DATA.completestr)));
            if step < 0
                id = id(id < DATA.commandctr);
                if ~isempty(id)
                    DATA.commandctr = id(end);
                else 
                    DATA.commandctr = 0;
                end
            else
                id = id(id > DATA.commandctr);
                if ~isempty(id)
                    DATA.commandctr = id(1);
                end
            end
        end
%        fprintf('Step Command %d\n',step);
        DATA.historyctr = length(DATA.oldcmds)+1;
        if DATA.commandctr == 0
            txt = DATA.oldcmds{end};
            DATA.commandctr = 1;
            DATA.historyctr = length(DATA.oldcmds);
        else
            txt = DATA.commands{DATA.commandctr};
            if DATA.commandctr < nlines
                if length(DATA.commandlines) >= DATA.commandctr
                    set(DATA.txtrec,'value',DATA.commandlines(DATA.commandctr));
                else
                    set(DATA.txtrec,'value',DATA.commandctr+1);
                end
            else
                set(DATA.txtrec,'value',nlines);
            end
        end
    elseif DATA.commandctr == 1 && (step < 0 || DATA.historyctr <= length(DATA.oldcmds))
        if DATA.historyctr == length(DATA.oldcmds) && step == 1
            DATA.historyctr = length(DATA.oldcmds)+1;
            txt = DATA.commands{1};
        elseif DATA.historyctr > 1
            if isempty(DATA.completestr) || DATA.historyctr == length(DATA.oldcmds)+1
                DATA.historyctr = DATA.historyctr+step;
            else
                id = find(strncmp(DATA.completestr,DATA.oldcmds,length(DATA.completestr)));
                if step < 0
                    id = id(id < DATA.historyctr);
                    if ~isempty(id)
                        DATA.historyctr = id(end);
                    end
                else
                    id = id(id > DATA.historyctr);
                    if ~isempty(id)
                        DATA.historyctr = id(1);
                    else
                        DATA.historyctr = 0;
                    end
                end
            end
            set(DATA.txtrec,'value',1);
            if DATA.historyctr == 0
                txt = 'Todays Commands';
                DATA.historyctr = length(DATA.oldcmds)+1;
                DATA.commandctr = 0;
            else
                txt = DATA.oldcmds{DATA.historyctr};
            end
        end
    else
        txt = src.Text;
    end
 %   txt = StripComments(txt);
    
function jTextKey(src, ev)    
    DATA = GetDataFromFig(src);
    ks =get(ev);
    newchar = 0;
    if ks.KeyCode == 38  %up arrow
        if ~isempty(DATA.completions)
            x = get(DATA.txtrec,'value');
            if x > 2
                set(DATA.txtrec,'value',x-1);
                src.Text = [DATA.completions{x-2} '='];
            end
        else
            [DATA, src.Text] = PrevCommand(DATA, src, -1);
        end
        src.setForeground(java.awt.Color(0.6,0,0))
        src.CaretPosition = length(src.Text);
    elseif ks.KeyCode == 10  %return
        if ~isempty(DATA.completions)
            DATA = ResetTextLst(DATA);
            src.setForeground(java.awt.Color(0,0,0));
            if isempty(strfind(src.Text,'='))  %User did not add to completion
                set(DATA.toplevel,'UserData',DATA);
            return;
            end
        else
        DATA = TextEntered(src, ev);
        end
    elseif ks.KeyCode == 39  %right arrow
    elseif ks.KeyCode == 40  %down arrow
        if ~isempty(DATA.completions)
            x = get(DATA.txtrec,'value');
            if x <= length(DATA.completions)
                set(DATA.txtrec,'value',x+1);
                src.Text = [DATA.completions{x} '='];
                src.CaretPosition = length(src.Text);
            end
        else
            [DATA, src.Text] = PrevCommand(DATA, src, 1);
        end
        src.setForeground(java.awt.Color(0.6,0,0))
        src.CaretPosition = length(src.Text);
    elseif ismember(ks.KeyCode,[8 127]) && strcmp(ks.ShiftDown,'on');  %shift delete
        src.Text = '';
    elseif ismember(ks.KeyCode,[8 127]) && strcmp(ks.ControlDown,'on');  %control-delete
        src.Text = regexprep(src.Text,'=.*','=');
        src.CaretPosition = length(src.Text);                        
    elseif isempty(ks.KeyChar)
    elseif ks.KeyChar == 9 %Tab
        a = deblank(src.Text);
        if isempty(strfind(a,'='))  %complete codes
            DATA = ShowCompletions(DATA,a);
        end
    elseif 0 && ks.KeyChar == ' ' %%now tab works, treat space normally
        a = deblank(src.Text);
        if isempty(strfind(a,'='))  %complete codes
            DATA = ShowCompletions(DATA,a);
        end
    else
        newchar = 1;
        if DATA.newchar == 0 %was loooking at history/completions
            txt = StripComments(src.Text);
            if ~strcmp(txt, src.Text)
                src.Text = txt;
                src.CaretPosition = length(src.Text);
                newchar = 0; 
            end
            src.setForeground(java.awt.Color(0,0,0));
        end
        DATA.completestr = '';
    end
    DATA.newchar = newchar; %something typed
    set(DATA.toplevel,'UserData',DATA);
        
function TextKey(src,ev)
    DATA = GetDataFromFig(src);
    ks.Key = get(ev,'KeyChar');
    if strcmp(ks.Key,'downarrow')
        if ~isempty(DATA.completions)
            x = get(DATA.txtrec,'value');
            if x <= length(DATA.completions)
                set(DATA.txtrec,'value',x+1);
                set(src,'string',[DATA.completions{x} '=']);
            end
        elseif DATA.commandctr < length(DATA.commands)
            DATA.commandctr = DATA.commandctr+1;
            set(src,'string',DATA.commands{DATA.commandctr})
            set(DATA.toplevel,'UserData',DATA);
        end
    elseif strcmp(ks.Key,'downarrow')
    elseif strcmp(ks.Key,'uparrow')
        if ~isempty(DATA.completions)
            x = get(DATA.txtrec,'value');
            if x > 3
                set(DATA.txtrec,'value',x-1);
                set(src,'string',[DATA.completions{x-2} '=']);
            end
        elseif DATA.commandctr > 1
            DATA.commandctr = DATA.commandctr-1;
            set(src,'string',DATA.commands{DATA.commandctr})
            set(DATA.toplevel,'UserData',DATA);
        end
    elseif strcmp(ks.Key,'space')
        a = deblank(get(src,'string'))
        if isempty(strfind(a,'='))  %complete codes
            DATA = ShowCompletions(DATA,a);
        end
    elseif strcmp(ks.Key,'tab')
        fprintf('Hit Tab');
    elseif strcmp(ks.Key,'control')
        setappdata(DATA.toplevel,'cntrl_is_down',1);
    else
        str = get(src,'string')
%        if str(end) ~= ' ' && ~isempty(DATA.completions)
%            ResetTextLst(DATA);
%        end
    end
    if strcmp(ks.Modifier,'control')
        setappdata(DATA.toplevel,'cntrl_is_down',1);
    else
        setappdata(DATA.toplevel,'cntrl_is_down',0);
    end
    
function DATA = ShowCompletions(DATA, a)
    DATA.completestr = a;
    xid = find(~strcmp('xx',{DATA.comcodes.code}));
        id = xid(find(strncmp(a,{DATA.comcodes(xid).code},length(a))));
        fprintf('Completing %s\n',a);
        str{1} = sprintf('%d Possible Completions for %s  (Return or Click here to cancel)',length(id),a);
        for j = 1:length(id)
            str{j+1} = [DATA.comcodes(id(j)).code '   (' DATA.comcodes(id(j)).label ')' ];
        end
% if completions not empty, there is already a set of completions up
% so don't change the recorded text
        if isempty(DATA.completions) %
            DATA.oldtxt = get(DATA.txtrec,'string');
        end
        DATA.completions = {DATA.comcodes(id).code};
        set(DATA.txtrec,'string',str,'foregroundcolor',[1 0 0]);
        if ~isempty(id)
            set(DATA.txtrec,'value',2);
        end
        set(DATA.toplevel,'UserData',DATA);

function DATA = ResetTextLst(DATA)       
    DATA.completions = {};
    set(DATA.txtrec,'string',DATA.oldtxt,'foregroundcolor',[0 0 0],'value',1);
    if nargout == 0
        set(DATA.toplevel,'UserData',DATA);
    end
    
function TextList(a,b)
    DATA = GetDataFromFig(a);
    line = get(a,'value');
    s = get(a,'string');
    if length(DATA.completions) >= line-1 && line > 1
        SetTextUI(DATA,[DATA.completions{line-1} '=']);
    elseif line ==1 && isfield(DATA,'oldtxt')
        ResetTextLst(DATA);
    else
        str = StripComments(s(line,:));
        SetTextUI(DATA, str);
    end

    
function s = StripComments(str)
    id = strfind(str,'(');
    if ~isempty(id)
        str = deblank(str(1:id(1)-1));
    end
    id = strfind(str,'?');
    if length(id) > 1
        str = deblank(str(2:id(2)-1));
    end
    s = regexprep(str,'\s\#.*','');

    
function SetTextUI(DATA, str)
if sum(strncmp(class(DATA.txtui),{'javahandle' 'jcontrol'},8))
    set(DATA.txtui,'Text', str);
else
    set(DATA.txtui,'string',str);        
end

    function yn = isjava(a)
        yn = 0;
        if isfield(a,'uipanel')
            yn = 1;
        end
        

function DATA = TextEntered(a,b)
    
    
    DATA = GetDataFromFig(a);
    if strncmp(class(a),'javahandle',8)
        txt = a.Text;
    else
        txt = get(a,'string');
        if get(gcf, 'currentcharacter') ~= 13 %return
            return;
        end
    end
    cntrl_is_down = getappdata(DATA.toplevel,'cntrl_is_down');
    if cntrl_is_down
        DATA = ShowCompletions(DATA,txt);
        return;
    elseif ~isempty(DATA.completions)
        DATA = ResetTextLst(DATA);
        if isempty(strfind(txt,'='))  %User did not add to completion
            set(DATA.toplevel,'UserData',DATA);
            return;
        end            
%        return;
    end
if isempty(txt)
return;
end

paused = PauseRead(DATA,1);
id = strfind(txt,'=');
if isstrprop(txt(1),'digit') || txt(1) == '-'
    txt = [DATA.lastcmd txt];
elseif length(id)
    DATA.lastcmd = txt(1:id(1));
end
    if DATA.currentstim > 1
        outprintf(DATA,'mo=%s\n',DATA.stimlabels{DATA.currentstim});
    end

    DATA = LogCommand(DATA, txt);
addline = 1;
str = [];
xstr = {};
if id
    code = txt(1:id(1)-1);
    id = strmatch(code,{DATA.comcodes.code},'exact');
    if length(id)
        str = DATA.comcodes(id(1)).label;
    else
    id =  strmatch(code,{DATA.strcodes.code},'exact');
    if length(id)
        str = DATA.strcodes(id(1)).label;
    end
    end
else
    [code,id] = FindCode(DATA,txt);
    if id
    str = DATA.comcodes(id).label;
    end
end

SetTextUI(DATA,'');
showbinoc = 0; % display value binoc returns

%Currently, binoc does NOT echo back everything verg sends.
%So only expect something if input produces a response
readargs = {};
if txt(end) ~= '='
    [DATA, a, b] = InterpretLine(DATA,txt,'fromgui');
    if a < 0
        addline = 0;
        showbinoc=2;
    end
elseif sum(strcmp(code,DATA.vergonlycodes))
    txt = [txt num2str(DATA.(code))];
else
    outprintf(DATA,'%s\n',txt);    
    readargs = {readargs{:} 'expect'};
    showbinoc = 1;
end
if ~isempty(regexp(txt,'=[A-z]+[+-]')) && ~isempty(code) %what is this for?
    outprintf(DATA,'%s=\n',code);
    showbinoc = 2;
    readargs = {readargs{:} 'expect'};
end
DATA = ReadFromBinoc(DATA,'from TextEntered ',readargs{:});
if showbinoc
%    code = txt(1:end-1);
    if sum(strcmp(code,{'uf' 'monkey'})) && isfield(DATA,'cwd')
       DATA = AddTextToGui(DATA,['cwd=' DATA.cwd]);
    end
    if strcmp(code,'cwd') && isfield(DATA,'cwd')
       DATA = AddTextToGui(DATA,['cwd=' DATA.cwd ': Set When set monkey=xxx']);
    end
    if sum(strcmp(code,{'rw'})) %show total reward
        if isfield(DATA.Trials,'good') && isfield(DATA.Trials,'rw')
            rwsum = sum([DATA.Trials([DATA.Trials.good] ==1).rw]);
            DATA = AddTextToGui(DATA,['totalreward=' sprintf('%.1f',rwsum)]);
        end
    end
    if strcmp(code,'op')
       txt = ['?' CodeText(DATA, 'optionflag')]; 
       str = 'Optionflag';
    elseif isfield(DATA.binoc{DATA.currentstim},code)
       id = strmatch(code,{DATA.comcodes.code},'exact');
       if length(id) == 1 && DATA.comcodes(id).type == 'C'
           txt = ['?' txt '?' DATA.binoc{DATA.currentstim}.(code)];
       elseif isempty(id)
           txt = ['?' txt '?' num2str(DATA.binoc{DATA.currentstim}.(code)') '(Unrecognized code)'];
       elseif ischar(DATA.binoc{DATA.currentstim}.(code)) %sometimes nmes->char  (!! ei=180lin)
           txt = ['?' txt '?''' DATA.binoc{DATA.currentstim}.(code) ''''];
       else
           if showbinoc ==2
               txt = [code '=' num2str(DATA.binoc{DATA.currentstim}.(code)')];
           elseif strncmp([code '='],DATA.lastline,length(code)+1)
               txt = ['?' DATA.lastline];
           else
               txt = ['?' txt '?' num2str(DATA.binoc{DATA.currentstim}.(code)')];
           end
       end
       if length(id)
        str = [DATA.comcodes(id(1)).label ',#' num2str(DATA.comcodes(id(1)).const-1)];
       end
    elseif strmatch(code,{'nr' 'nt' 'n2' 'n3' 'et' 'e2' 'e3' })
       txt = ['?' CodeText(DATA, code)];       
       str = 'Nstim';
    elseif isfield(DATA.binocstr,code)
        txt = ['?' txt '?''' DATA.binocstr.(code) ''''];
       id =  strmatch(code,{DATA.strcodes.code},'exact');
       if length(id)
           str = DATA.strcodes(id(1)).label;
       end
    else
        id = find(strncmp(code,{DATA.comcodes.code},length(code))); %find possible matches
        if ~isempty(id)
            str = 'Poossible Completions';
            for j = 1:length(id)
                xstr{j} = sprintf('%s %s\n',DATA.comcodes(id(j)).code,DATA.comcodes(id(j)).label);
            end
        else
            txt = [txt ' ' DATA.lastline];
        end                
    end
end


if addline %show this instruction in the history window. Cancelled commands dont
    a =  get(DATA.txtrec,'string');
    txtpos =  get(DATA.txtrec,'value');
    n = size(a,1);
    DATA.commandlines(length(DATA.commands)) = n+1;
    if strcmp(code,'px')
        pixdeg = atan(DATA.binoc{1}.px/DATA.binoc{1}.vd) * 180/pi;
        txt = [txt sprintf(' (%.5f deg)',pixdeg)];
    end
    if length(str)
        txt = [txt '(' str ')  ' datestr(now,'HH:MM')];
    end
    if iscellstr(a)
        a{n+1} = txt;
    else
        a(n+1,1:length(txt)) = txt;
    end
    for j = 1:length(xstr)
        a(n+1+j,1:length(xstr{j})) = xstr{j};
    end
    set(DATA.txtrec,'string',a);
    x = get(DATA.txtrec,'value');
%used only to do this if x > size(a,1) but need
%to set value to length a to make added text visible
%so added the elseif...
    if x > size(a,1)
        set(DATA.txtrec,'value',size(a,1));
    elseif x > 1
        set(DATA.txtrec,'value',size(a,1));        
    end
    set(DATA.txtrec,'listboxtop',n+1);
    DATA = LogCommand(DATA, txt, 'norec');
end
set(DATA.toplevel,'UserData',DATA);
SetGui(DATA);
if paused == 0 %wasn't paused at start
    paused = PauseRead(DATA,0);
end


function ts = FindSessionStart(DATA)
    cmdfile = ['/local/' DATA.binoc{1}.monkey '/binoccmdhistory'];
    txt = scanlines(cmdfile);
    s = datestr(now,'dd-mmm-yyyy');
    id = find(CellToMat(strfind(txt,s)));
    ts = 0;
    j = 1;
    while ts == 0 && j < length(id)
        d = strfind(txt{id(j)},s); 
        if sum(strncmp(txt{id(j)},{'# Run Hit'},6))
            ts = datenum(txt{id(j)}(d:end));
        end
        j = j+1;
    end
    
     txt = scanlines(['/local/' DATA.binoc{1}.monkey '/pen' num2str(DATA.binoc{1}.Pn) '.log']);
     id = find(CellToMat(strfind(txt,'Entered Brain')));
     id = find(CellToMat(strfind(txt,'Head Restrained')));
              

function DATA = AddTextToGui(DATA, txt, varargin)
    if ~isfield(DATA,'txtrec') || ~ishandle(DATA.txtrec)
        return;
    end
a =  get(DATA.txtrec,'string');
n = size(a,1);
DATA = LogCommand(DATA, txt, varargin{:});
txt  = [txt ' ' datestr(now,'HH:MM')];
a(n+1,1:length(txt)) = txt;
set(DATA.txtrec,'string',a);
set(DATA.txtrec,'listboxtop',n+1);

        

function ClearTaggedChecks(it, tags)

    
    c = get(it,'children');
    menutags = get(c,'tag');
    if isempty(tags)
        for j = 1:length(menutags)
            if ~isempty(menutags{j})
                set(c(j),'Checked','off');
            end
        end
    else
        for j = 1:length(tags)
            k = strmatch(tags{j},menutags,'exact');
            set(c(k),'Checked','off');
        end
    end
    
function ChoosePsych(a,b, mode)
    DATA = GetDataFromFig(a);
    onoff = {'off' 'on'};
    F = GetFigure(a);
    Expts = getappdata(F,'Expts');
    if strncmp(mode,'Expt',4)
        e = sscanf(mode,'Expt%d');
        if length(DATA.plotexpts) < e
            DATA.plotexpts(e) = 1;
        else
            DATA.plotexpts(e) = ~DATA.plotexpts(e);
        end
        if strcmp(DATA.psych.blockmode,'OneOnly')
            DATA.plotexpts = zeros(size(DATA.plotexpts));
            DATA.plotexpts(e) = 1;
            c = get(get(a,'parent'),'children');
            for j = 1:length(c)
                if ~strcmp(get(c(j),'tag'),'OneOnly')
                    set(c(j),'checked','off');
                end
            end
        else
            ClearTaggedChecks(get(a,'parent'),{});
        end
        if strmatch(DATA.psych.blockmode,{'Current' 'OneOnly'})
            it = findobj(get(a,'parent'),'Tag',DATA.psych.blockmode);
            set(it,'Checked','on');
        else
            DATA.psych.blockmode = 'Select';
        end
        PlotPsych(DATA, Expts);

        set(a,'checked',onoff{DATA.plotexpts(e)+1});
    elseif strcmp(mode,'exptrelist')
        PsychMenu(DATA);
    elseif strcmp(mode,'exptsummary')
        PlotExptsSummary(DATA.Expts);
    elseif strmatch(mode,'Current')
        if strcmp(DATA.psych.blockmode,'Current')
            DATA.psych.blockmode = 'Select';
            set(a,'checked','off');
        else
            DATA.psych.blockmode = mode;
            set(a,'checked','on');
        end
        ClearTaggedChecks(get(a,'parent'),{});
        PlotPsych(DATA);
    elseif strmatch(mode,{'OnlyCurrent','All' 'None' 'OneOnly'})
        DATA.plotexpts = zeros(1,length(DATA.expts));
        c = get(get(a,'parent'),'children');
        set(c,'checked','off');
        if ~strcmp(DATA.psych.blockmode,mode)
            set(a,'checked','on');
            DATA.psych.blockmode = mode;
        else
            set(a,'checked','on');
            DATA.psych.blockmode = 'Select';
        end        
        if strcmp(mode,'All')
            DATA.plotexpts(1:end) = 1;
            strs = get(c,'label');
            id = find(strncmp('Expt',strs,4));
            set(c(id),'checked','on');
        end
        PlotPsych(DATA);
    elseif strmatch(mode,'Pause')
        DATA.psych.show = ~DATA.psych.show;
        set(a,'Checked',onoff{DATA.psych.show+1});
        PlotPsych(DATA);
    elseif sum(strcmp(mode,{'collapse2' 'collapse3'}))
        j = find(strcmp(mode,{'collapse2' 'collapse3'}))+1;
        DATA.psych.collapse(j) = ~DATA.psych.collapse(j);
        set(a,'Checked',onoff{DATA.psych.collapse(j)+1});
        PlotPsych(DATA);
    elseif sum(strcmp(mode,{'crosshairs' 'trialresult' 'showblocks'}))
        DATA.psych.(mode) = ~DATA.psych.(mode);
        set(a,'Checked',onoff{DATA.psych.(mode)+1});
        PlotPsych(DATA);
    elseif strmatch(mode,'choosepsychfile')
        [name, path] = uigetfile(['/local/data/psych/' DATA.binoc{1}.monkey '/*.*'],'select psych data file');
        if name
            Expts = ReadPsychFile([path '/' name]);
            setappdata(gcf,'Expts',Expts);
            PsychMenu(DATA,Expts);
            GetFigure('Expts Summary');
            PlotExptsSummary(Expts);
        end
    elseif strmatch(mode,'readpsychfile')
        logfile = sprintf('/local/%s/logs/%s%s',DATA.binoc{1}.monkey,DATA.binoc{1}.monkey,datestr(now,'ddmmmyyyy')); 
        if exist(DATA.binoc{1}.psychfile,'file')
            Expts = ReadPsychFile(DATA.binoc{1}.psychfile);
        elseif exist(logfile,'file')
            Expts = ReadPsychFile(logfile,'useallexpts','nmin',3);
        else            
            Expts = {};
        end
        if ~isempty(Expts)
            setappdata(gcf,'Expts',Expts);
            PsychMenu(DATA,Expts);
        end
    elseif strcmp(mode,'savetrials');
        [outname, pathname] = uiputfile(['/local/' DATA.binoc{1}.monkey '/PsychDat.mat']);
        if outname
            outname = [pathname '/' outname];
            Data = rmfields(DATA,{'timerobj' 'txtui'})'
            save(outname,'Data');
            fprintf('Trials saved to %s\n',outname);
        end
    end
    set(DATA.toplevel,'UserData',DATA);



function DATA = SetFigure(tag, DATA)

    [a,isnew] = GetFigure(tag,'parent',DATA.toplevel);
    onoff = {'off' 'on'};
    if isnew
        DATA.figs.(tag) = a;
        if strcmp(tag,'VergPsych')
            hm = uimenu(a, 'Label','Expts','Tag','ExptMenu');
            PsychMenu(DATA);
            hm = uimenu(a, 'Label','Options','Tag','PsychOptions');
            sm = uimenu(hm,'Label','Crosshairs','callback', {@ChoosePsych, 'crosshairs'},...
                'checked',onoff{DATA.psych.crosshairs+1});
            sm = uimenu(hm,'Label','Collapse Expt2','callback', {@ChoosePsych, 'collapse2'},...
                'checked',onoff{DATA.psych.collapse(2)+1});
            sm = uimenu(hm,'Label','Collapse Expt3','callback', {@ChoosePsych, 'collapse3'},...
                'checked',onoff{DATA.psych.collapse(3)+1});
            sm = uimenu(hm,'Label','Just Show Trial outcomes','callback', {@ChoosePsych, 'trialresult'},...
                'checked',onoff{DATA.psych.trialresult+1});
            sm = uimenu(hm,'Label','Separate by Block','callback', {@ChoosePsych, 'showblocks'},...
                'checked',onoff{DATA.psych.showblocks+1});
            sm = uimenu(hm,'Label','Save Trial Data','callback', {@ChoosePsych, 'savetrials'});
            sm = uimenu(hm,'Label','Read Todays PsychFile','callback', {@ChoosePsych, 'readpsychfile'});
            sm = uimenu(hm,'Label','Read Previous PsychFile','callback', {@ChoosePsych, 'choosepsychfile'});
            sm = uimenu(hm,'Label','Plot Expt Summary','callback', {@ChoosePsych, 'exptsummary'});
            sm = uimenu(hm,'Label','Refresh Expt List','callback', {@ChoosePsych, 'exptrelist'});
            set(a,'UserData',DATA.toplevel);
            set(a,'DefaultUIControlFontSize',DATA.font.FontSize);
            set(a,'DefaultUIControlFontName',DATA.font.FontName);

        end
        set(DATA.toplevel,'UserData',DATA);
    end

function PsychMenu(DATA, varargin)    
    
    Expts = DATA.Expts;
    j = 1; 
    while j <= length(varargin)
        if iscell(varargin{j})
            Expts = varargin{j};
        end
        j = j+1;
    end
    if ~isfield(DATA,'figs') || ~isfigure(DATA.figs.VergPsych)
        return;
    end
    hm = findobj(DATA.figs.VergPsych,'tag','ExptMenu');
    c = get(hm,'children');
    delete(c);
    for j = 1:length(Expts)
        if isfield(Expts{j},'Trials')
        nt = length(Expts{j}.Trials);
        sm = uimenu(hm,'Label', sprintf('Expt%d %s %d',j,Expt2Name(Expts{j}),nt),'CallBack', {@ChoosePsych, sprintf('Expt%d',j)});
        if j < length(DATA.plotexpts) && DATA.plotexpts(j)
            set(sm,'Checked','on');
        end
        end
    end
    sm = uimenu(hm,'Label', 'Current','Callback',{@ChoosePsych, 'Current'});
    if strcmp(DATA.psych.blockmode,'Current')
        set(sm,'Checked','On');
    end
    uimenu(hm,'Label', 'Only Current','Callback',{@ChoosePsych, 'OnlyCurrent'},'tag','OnlyCurrent');
    uimenu(hm,'Label', 'All','Callback',{@ChoosePsych, 'All'},'tag','All');
    uimenu(hm,'Label', 'One','Callback',{@ChoosePsych, 'OneOnly'},'tag','OneOnly');
    uimenu(hm,'Label', 'None','Callback',{@ChoosePsych, 'None'},'tag','None');
    
function DATA = CheckExpts(DATA)

    for j = 1:length(DATA.Expts)
        if j < length(DATA.Expts) && ~isfield(DATA.Expts{j},'last')
            DATA.Expts{j}.last = DATA.Expts{j+1}.first -1;
        end
        DATA.Expts{j}.Header.exptno = j;
        DATA.Expts{j}.Header.rc = 0; %Cant' do rc online at the moment
    end
    
function DATA = PlotPsych(DATA, Expts)

    forcewin = 0;

    if nargin ==1
        Expts = {};
    elseif ischar(Expts)
        vararg = Expts;
        Expts = {};
        if strcmp(vararg,'force');
            forcewin =1;
        end
    end
    
    DATA.err = '';
    if (isempty(DATA.Expts) || isempty(DATA.Trials)) && isempty(Expts)
        DATA.err = sprintf('No Trials/Expts');
%        return;
    end
    if strmatch(DATA.psych.blockmode,'None')
        DATA.err = sprintf('Plotting Off - See ''Expts'' Menu');
        return;
    end
    if isempty(Expts)
        DATA = CheckExpts(DATA);
        Expts = DATA.Expts;
        

    Expt = [];
    e = length(DATA.Expts);
    if e == 0 && forcewin ==0
        return;
    end
    
    allid = [];
    if strcmp(DATA.psych.blockmode,'All') || sum(DATA.plotexpts)
        if strmatch(DATA.psych.blockmode,'All')
            expts = 1:e-1;
        else
            expts = find(DATA.plotexpts);
        end
        %can only combine expts if have same et,e2 types.  Check each
        %against last in list
        if sum(strcmp(DATA.psych.blockmode,{'Current' 'All'}))
            Expt = Expts{e};
        else
            Expt = Expts{expts(end)};
        end
        
        for j = expts
            if isfield(Expts{j},'last')
                last = Expts{j}.last;
            else
                last = length(DATA.Trials);
            end                
            if strcmp(Expts{j}.Stimvals.et,Expt.Stimvals.et) && ...
               strcmp(Expts{j}.Stimvals.e2,Expt.Stimvals.e2)
                allid= [allid Expts{j}.first:last];
            end
        end
    elseif sum(e) > 0
        Expt = Expts{e};
    end
    if forcewin
        DATA = SetFigure('VergPsych', DATA);
    end
    if isempty(Expt)
        return;
    end

    %Always update trial list for current expt
    id = DATA.Expts{e}.first:length(DATA.Trials);
        DATA.Expts{e}.Trials = DATA.Trials(id);
    if Expt.first == DATA.Expts{e}.first  %using current
        if strmatch(DATA.psych.blockmode,{'All' 'Current' 'OnlyCurrent'})
         allid = [allid id];
        end
    end
    id = unique(allid);
    if length(id) < 2
        DATA.err = sprintf('Too Few Trials(%d)',length(id));
        return;
    end
    Expt.Trials = DATA.Trials(id);
    Expt.Header.rc = 0;
    Expt.Header.expname  = 'Online';
    if isfield(Expt.Trials,'RespDir')
        Expt.Header.psych  = 1;
    else
        Expt.Header.psych  = 0;
    end
    
    f = {'m2' 'or' 'm3' 'n2'};
    for j = 1:length(f)
        if isfield(DATA.binoc{1},f{j})
            Expt.Stimvals.(f{j}) = DATA.binoc{1}.(f{j});
        end
    end
    else
        expts = find(DATA.plotexpts);
        Expt = CombineExpts(Expts(expts));
    end    
    
    if DATA.psych.show
        for j = 1:length(Expt.Trials)
            if ~isempty(Expt.Trials(j).Trial)
                good(j) = 1;
            else
                good(j) = 0;
            end
        end
        Expt.Trials = Expt.Trials(find(good));
                
        DATA = SetFigure('VergPsych', DATA);
        hold off;
        eargs = {'collapse' DATA.psych.collapse};
        if DATA.psych.trialresult
            plot([Expt.Trials.Start],[Expt.Trials.good],'o');
            datetick('x','hh:mm');
            axis('tight');
        elseif DATA.psych.showblocks
            ExptPsych(Expts(expts),'labelblock','nmin',2, 'mintrials', 2, 'shown', eargs{:});
        else
        np = sum(abs([Expt.Trials.RespDir]) ==1); %psych trials
        Expt = FillTrials(Expt,Expt.Stimvals.et);
        Expt = FillTrials(Expt,Expt.Stimvals.e2);
        if isfield(Expt.Trials,'psyv')
            npsy = length(unique([Expt.Trials.psyv]));
            if npsy > 2
                eargs = {eargs{:} 'type' 'psyv'};
            end
        end
        if isfield(DATA.matexpres,'types') && length(DATA.matexpres.types) > 1
            eargs = {eargs{:} 'type2' DATA.matexpres.types{2}};
        elseif ~strcmp(Expt.Stimvals.e2,'e0')
            eargs = {eargs{:} 'type2' Expt.Stimvals.e2};
        end
        if np > 1
            if strcmp(Expt.Stimvals.e2,'od') && isfield(Expt.Stimvals,'or')
                for j = 1:length(DATA.expvals{2})
                    do = DATA.expvals{2}(j);
                    lo = Expt.Stimvals.or + do/2;
                    ro = Expt.Stimvals.or - do/2;
                    legendlabels{j} = sprintf('R%.0dL%.0f',ro,lo);
                end
                eargs = {eargs{:} 'legendlabels' legendlabels};
            end
            try
                [a,b] = ExptPsych(Expt,'nmin',1,'mintrials',2,'shown',eargs{:});
            catch ME
                DATA.err = sprintf('Error In ExptPsych %s',ME.message);
                fprintf('%s\n',DATA.err);
            end
            id = find(strcmp(Expt.Stimvals.et,{DATA.comcodes.code}));
            if length(id) == 1
                set(get(gca,'xlabel'),'string',DATA.comcodes(id).label);
            end
            if DATA.psych.crosshairs
                holdstate = ishold;
                hold on;
                plot([0 0], get(gca,'ylim'),'k:');
                plot(get(gca,'xlim'),[0.5 0.5],'k:');
                if holdstate == 0
                    hold off;
                end
            end
        end
        end
    end
    DATA.Expt = Expt;
    
    
function DATA = CheckTimer(DATA)
    
    global rbusy;
    
    if isfield(DATA,'timerobj') && isobject(DATA.timerobj)
        on = get(DATA.timerobj,'Running');
    end

    DATA.pipeerror = 0;
    c = get(DATA.toplevel,'color');
    if rbusy > 0
        set(DATA.toplevel,'color','r');
        DATA.pipeerror = 1;
    end
    if strcmp(get(DATA.timerobj,'running'),'off')
        set(DATA.toplevel,'color',[0 0 0.5]);
        DATA.pipeerror = 2;
    end
    if DATA.pipeerror == 0 && sum(c == DATA.windowcolor) < 3
        set(DATA.toplevel,'color',DATA.windowcolor);        
    end
