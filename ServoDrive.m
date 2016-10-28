function DATA = ServoDrive(varargin)
%DATA = ServoDrive(...
%Matlab controller for Servo controlled microdrive, using serial port
%ServoDrive('ttyname', device)
%sets the name of the serial port connected to use.
%This can also be set by creating a file '/local/servomotor.setup' with a line
%    serialport=device
%
%
%ServoDrive(...,'position',x) uses x for the figure position.
%
%ServoDrive(...,'callback',fcn) gives a callback handle. This will
%be called at the end of each movement with the new position
%
%DATA = ServoDrive('readpos') returns a data structure with the state
%   including  DATA.position (current depth in microns)
%              DATA.alldepths
%              DATA.alltimes    depths set since starting, and their times

figpos = [1400 400 350 350];
DATA.figpos = [];
tag = 'Servo Controller';
DATA.verbose = 0;
verbose = [];
startdepth = NaN;


%
% Command SR1 -> SR20  changes period with which servo is polled (10KHz ->
% 500Hz). But does not seem likely to help with noise
%
%Seems clear the encoder on the motor just measures displacement (its an OSE), and
%the chip intergrates this. So disconnecting cable to motor, then moving
%the drive and reconnecting proudces NO CHANGE in sensed position.
%Conversely powering off the chip causes the position to be lost/reset. 

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'callback',5)
        j = j+1;
        DATA.callback = varargin{j};
        if length(DATA.callback) > 1 && isfigure(DATA.callback{2})
            DATA.callingfigure = DATA.callback{2};
        end
    elseif strncmpi(varargin{j},'depths',5)
        j = j+1;
        DATA.alldepths = varargin{j};
        j = j+1;
        DATA.alltimes = varargin{j};
        DATA = CheckDepthData(DATA);
    elseif strncmpi(varargin{j},'position',5)
        j = j+1;
        DATA.figpos = varargin{j};
    elseif strncmpi(varargin{j},'quiet',5)
        DATA.verbose = 0;
        verbose = 0;
    elseif strncmpi(varargin{j},'startdepth',8)
        j = j+1;
        startdepth = varargin{j};
    elseif strncmpi(varargin{j},'stepsize',5)
        j = j+1;
        DATA.stepsize = varargin{j};
    elseif strncmpi(varargin{j},'ttyname',3)
        j = j+1;
        DATA.ttyname = varargin{j};
    elseif strncmpi(varargin{j},'verbose',5)
        DATA.verbose = 1;
        verbose = 1;
    end
    j = j+1;
end

[F, isnew] = GetFigure(tag);
if isnew
    DATA.toplevel = F;
    DATA = SetDefaults(DATA);
    set(F,'position',DATA.figpos);
    DATA = BuildServoWindow(DATA);
    DATA = OpenServoPort(DATA);
    DATA = OpenDiskLog(DATA);
    DATA = SetStartDepth(DATA, startdepth);
    set(DATA.toplevel,'UserData',DATA);
else 
    DATA = get(F,'UserData');
end
if ~isempty(verbose) %forced on command line
    DATA.verbose = verbose;
end

%Things to do after startup
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'readpos',5)
        d = GetCurrentPosition(DATA, 'show');
        DATA.alldepths = [DATA.alldepths d];
        DATA.alltimes = [DATA.alltimes now];
        DATA.position = d;
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(varargin{j},'label',5)
        set(DATA.toplevel,'Name',varargin{j+1});
    elseif strncmpi(varargin{j},'close',5)
        CloseServoPort(DATA);

    end
    j = j+1;
end

function MoveMicroDrive(a,b,fcn)
DATA = GetDataFromFig(a);

EnableMoveButtons(DATA,'off');
it = findobj(DATA.toplevel,'tag','StepSize');

if strcmp(fcn,'moveto')
    newpos = str2num(get(DATA.setdepth,'string'));
    v = rpm2speed(DATA);
%For large moves at 100uM/sec or faster, check with the user    
    if abs(newpos - DATA.position) > 1000 && v > 0.09
        str = sprintf('Do you want speed of %.3fuM/s, or automatic (slow)?',v.*1000);
        yn = questdlg(str,'Large Fast Move', 'Yes', 'Auto','Yes');
        if ~strcmp(yn,'Yes')
            DATA.lastspeed = DATA.motorspeed;
            DATA.motorspeed = 0; %force automatic speed setting
            SetData(DATA);
        end
    end
    DATA = SetNewPosition(DATA, newpos);
elseif strcmp(fcn,'readposition')
    d = GetCurrentPosition(DATA, 'show');
    DATA.position = d;
elseif strcmp(fcn,'set')
    newpos = str2num(get(DATA.setdepth,'string'));
    SetHomePosition(DATA, newpos);
else
    str = get(it,'string');
    j= get(it,'value');
    step = sscanf(str(j,:),'%d');
    DATA.stepsize = step;
    if  strcmp(fcn,'Up')
        step = -step;
    end
    if DATA.verbose > 1
        fprintf('Moving %d\n',step);
    end
    DATA = ChangePosition(DATA, step);
end
EnableMoveButtons(DATA,'on');
set(DATA.toplevel,'UserData',DATA);


function DATA = CheckDepthData(DATA)
if length(DATA.alltimes) > length(DATA.alldepths)
    DATA.alltimes = DATA.alltimes(1:length(DATA.alldepths));
end

function DATA = SetStartDepth(DATA, startdepth)

% set initial position, but check that drive has not already been used;
if isnan(startdepth)
%not setting startdepth. So if alldepths exists, this should be a continuation
%If depth does not match last recorded, Motor might have been turned off.
%Check whether to restore previous value

    yn = 'No';
    if isfield(DATA,'alldepths') && ~isempty(DATA.alldepths)
        startdepth = DATA.alldepths(end);
        x = GetCurrentPosition(DATA);
        if abs(startdepth-DATA.position) > 10
            msg = sprintf('MicroDrive was %.0fuM At %s. Now Reads %.0f\nPerhaps it was turned off.\n Reset (no movement) to  %.0f?',DATA.alldepths(end),datestr(DATA.alltimes(end)),x,startdepth);
            yn = questdlg(msg,'ServoDrive Message','Yes','No','Yes');
        end
    end
    if strcmp(yn,'No')
        return;
    end    
elseif length(DATA.alldepths) && max(DATA.alltimes)-now < 2/24
    x = GetCurrentPosition(DATA);
    msg = sprintf('MicroDrive last set to %.0fuM on %s. Now Reads %.0f. Do you you want to set %.0f (without moving)?',DATA.alldepths(end),datestr(DATA.alltimes(end)),x,startdepth);
    yn = questdlg(msg,'ServoDrive Message','Yes',sprintf('No Leave Reading at %.3f',x),'Yes');
    if strncmp(yn,'No',2)
        return;
    end    
end

newpos = startdepth;
SetHomePosition(DATA, newpos);

function DATA = BuildServoWindow(DATA)
F = DATA.toplevel;
nr = 7;
nc=6;

set(F,'menubar','none','numbertitle','off');
if exist('/local/Images/icons/DownArrow.mat')
    im = load('/local/Images/icons/DownArrow.mat');
else
    fprintf('Cant Read /local/Images/icons/DownArrow.mat. Using simple arrrow\n');
    im.x = zeros(75,75);
    [x ,y] = meshgrid([1:75]-37,1:75);
    ii  = find(y/1.5 + abs(x) < 50 & y > 15);
    im.x(ii) = 1;
end
im.x = im.x-min(im.x(:));
im.x = im.x./max(im.x(:));
[ii,ij] = ind2sub(size(im.x),find(im.x < 0.04));

im.x = 1-im.x;
upx = flipud(im.x);
bcolor = get(gcf,'color');
im.cdata = repmat(upx,[1 1 3]);
im.cdata(find(im.cdata > 0.96)) = mean(bcolor);
bp = [0.1 0.7 0.3 0.2];
uicontrol(gcf,'style','pushbutton','string','Up', ...
     'foregroundcolor','w', ...
   'Callback', {@MoveMicroDrive, 'Up'}, 'Tag','UpArrow',...
   'backgroundcolor',bcolor,...
        'units', 'norm', 'position',bp,'value',1,'cdata',im.cdata);


bp = [0.05 0.4 0.3 0.15];
h = uicontrol(gcf,'style','popupmenu','string','10uM|20uM|25uM|50uM|100uM|200uM|Set Custom', ...
        'units', 'norm', 'position',bp,'value',1,...
        'Tag','StepSize','fontsize',18,'fontweight','bold','callback',@SetStepSize);
if isfield(DATA,'stepsize') && DATA.stepsize > 0
    [a,b] = min(abs(DATA.stepsize-[10 20 25 50 100 200]));
    if a < 2
        set(h,'value',b');
    else
        DATA.customstep = DATA.stepsize;
        s = get(h,'string');
        val = size(s,1);
        str = sprintf('%d (Custom)',DATA.customstep);
        DATA.step = DATA.customstep;
        s(val+1,1:length(str)) = str;
        set(h,'string',s,'value',val+1);
    end
end
    bp(2) = 0.2;
    bp(1) = bp(1)+bp(3)+0.05;
bp(3) = 0.15;
uicontrol(gcf,'style','pushbutton','string','Set', ...
   'Callback', {@MoveMicroDrive, 'set'}, 'Tag','SetPosition',...
   'fontsize',18,'fontweight','bold',...
        'units', 'norm', 'position',bp,'value',1);
bp(1) = bp(1)+bp(3);
bp(3) = 0.2;
DATA.setdepth = uicontrol(gcf,'style','edit','string','0', ...
   'fontsize',18,'fontweight','bold',...
   'Tag', 'ManualSet',...
        'units', 'norm', 'position',bp,'value',1);
bp(1) = bp(1)+bp(3);
bp(3) = 0.2;
uicontrol(gcf,'style','pushbutton','string','Move', ...
   'Callback', {@MoveMicroDrive, 'moveto'}, 'Tag','MoveButton',...
   'fontsize',18,'fontweight','bold',...
        'units', 'norm', 'position',bp,'value',1);
% 
% bp = [0.5 0.2 0.45 0.15];
%     uicontrol(gcf,'style','pushbutton','string','Reopen', ...
%    'Callback', {@MoveMicroDrive, 'reopen'}, ...
%    'fontsize',18,'fontweight','bold',...
%         'units', 'norm', 'position',bp,'value',1);

bp = [0.1 0.21 0.3 0.2];
im.cdata = repmat(im.x,[1 1 3]);
im.cdata(find(im.cdata > 0.96)) = mean(bcolor);
h = uicontrol(gcf,'style','pushbutton','string','Down', ...
   'Callback', {@MoveMicroDrive, 'down'}, 'Tag','DownArrow',...
        'units', 'norm', 'position',bp,'value',1,...
        'foregroundcolor','w', 'backgroundcolor',bcolor,'cdata',im.cdata);
    bp(1) = bp(1)+bp(3);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.1;

bp = [0.05 0.1 0.9 0.15];
h = uicontrol(gcf,'style','check','string','Emergency Stop', ...
        'units', 'norm', 'position',bp,'value',0,...
        'Tag','EmergencyStop','fontsize',18,'fontweight','bold','callback',@StopMotor);

bp = [0.3 0.01 0.4 0.1];
DATA.depthlabel = uicontrol(gcf,'style','text','string','0 uM', ...
        'units', 'norm', 'position',bp,'value',1,'backgroundcolor',bcolor,'foregroundcolor',[0.8 0.3 0.2],...
        'fontsize',24,'fontweight','bold');
mn = uimenu(F,'Label','Close','callback',@CloseServoPort);
mn = uimenu(F,'Label','Reopen','callback',@OpenServoPort);
mn = uimenu(F,'Label','Stop','callback',@StopMotor,'Tag','StopButton');
mn = uimenu(F,'Label','Options');
sm = uimenu(mn,'Label','Read Position','callback',{@MoveMicroDrive, 'readposition'},'accelerator','R');
sm = uimenu(mn,'Label','Plots');
uimenu(sm,'Label','Full History','tag','FullHistory','checked','on','callback',@SetMotorPlot);
uimenu(sm,'Label','Last Move','tag','LastMove','callback',@SetMotorPlot);
uimenu(sm,'Label','Set Time range','tag','Custom','callback',@SetMotorPlot);
uimenu(sm,'Label','Speed','tag','MoveSpeed','callback',@SetMotorPlot);
uimenu(sm,'Label','None','tag','None','callback',@SetMotorPlot);
sm = uimenu(mn,'Label','Speed');
uimenu(sm,'Label','Normal (1 mm/s)','tag','Normal','Checked','on','callback',@SetMotorSpeed);
uimenu(sm,'Label','Fast (5 mm/s)','tag','Fast','callback',@SetMotorSpeed);
uimenu(sm,'Label','Slow (0.5 mm/s)','tag','Slow','callback',@SetMotorSpeed);
uimenu(sm,'Label','Very Slow (0.1 mm/s)','tag','VSlow','callback',@SetMotorSpeed);
uimenu(sm,'Label','Retracting (0.025 mm/s)','tag','VVSlow','callback',@SetMotorSpeed);
uimenu(sm,'Label','Automatic','tag','Auto','callback',@SetMotorSpeed);
uimenu(sm,'Label','Custom','tag','Custom','callback',@SetMotorSpeed);
sm = uimenu(mn,'Label','Verbose','callback',@SetOption,'tag','verbose');
if DATA.verbose
    set(sm,'checked','on');
end
sm = uimenu(mn,'Label','Show Buttons','callback',{@MenuOptions,'showbuttons'});
if DATA.verbose
    set(sm,'checked','on');
end

set(gca,'position',[0.4 0.4 0.6 0.6],'xtick',[],'ytick',[]);
set(DATA.toplevel,'UserData',DATA);

function EnableMoveButtons(DATA,state)

tags = {'DownArrow', 'UpArrow' 'MoveButton'};
for j = 1:length(tags)
    it = findobj(DATA.toplevel,'tag',tags{j});
    if ~isempty(it)
        set(it,'visible',state);
    end
end
it = findobj(DATA.toplevel,'tag','StopButton');
set(it,'enable','on');

function MenuOptions(a,b,fcn)
   DATA = GetDataFromFig(a);
if strcmp(fcn,'showbuttons')
    EnableMoveButtons(DATA,'on');
end



function SetOption(a,b)
   DATA = GetDataFromFig(a);
   tag = get(a,'tag');
   DATA.(tag) = ~DATA.(tag);
   if DATA.(tag) 
       set(a,'checked','on');
   else
       set(a,'checked','off');
   end
   set(DATA.toplevel,'UserData',DATA);
   
   
function DATA = OpenDiskLog(DATA)
%If the history file has positions from the last 12 hours
%read them into the history that will be displayed/checked
if exist(DATA.logfile,'file')
    X = load(DATA.logfile);
    logdate = max(X.alltimes);
    a = datevec(logdate);
    b = datevec(now);
    if now - logdate > 0.5  || a(3) < b(3)
        BackupFile(DATA.logfile);
    else
        X = load(DATA.logfile);
        DATA.alldepths = X.alldepths;
        DATA.alltimes = X.alltimes;
        DATA = CheckDepthData(DATA);
    end
end

function DATA = SaveDiskLog(DATA)

X.alltimes = DATA.alltimes;
X.alldepths = DATA.alldepths;
save(DATA.logfile,'-struct','X');


function SetMotorPlot(a,b)

   DATA = GetDataFromFig(a);
   tag = get(a,'tag');
   if strcmp(tag,'Custom')
       if DATA.timerange > 0
           defaultanswer{1} = num2str(DATA.timerange(1));
       else
           defaultanswer{1} = '10';
       end
       str = inputdlg('Plot How many minutes back?','Custom',1,defaultanswer);
       if isempty(str)
           return;
       end
       DATA.timerange = str2num(str{1});
   else
       DATA.plottype = tag;
       if sum(strncmp(tag,{'FullHistory' 'LastMove'},6))
           DATA.timerange = 0;
       end
       SetMenuCheck(a,'exclusive');
   end
   DATA.firstsample = 0;
   PlotDepths(DATA);
   set(DATA.toplevel,'UserData',DATA);

function SetMotorSpeed(a,b)
DATA = GetDataFromFig(a);
tag = get(a,'tag');

if strcmp(tag,'Normal')
    DATA.motorspeed = 1.0;
elseif strcmp(tag,'Slow')
    DATA.motorspeed = 0.5;
elseif strcmp(tag,'Fast')
    DATA.motorspeed = 5;
elseif strcmp(tag,'VSlow')
    DATA.motorspeed = 0.1;
elseif strcmp(tag,'VVSlow')
    DATA.motorspeed = 0.025;
elseif strcmp(tag,'Auto')
    DATA.motorspeed = 0;
elseif strcmp(tag,'Custom')
    defaultanswer{1} = num2str(DATA.customspeed.*1000);
    str = inputdlg('Step Size in uM/sec','Custom',1,defaultanswer);
    if isempty(str)
       return;
    end
    step = str2num(str{1});
    str = sprintf('Custom (%d uM/s)',step);
    set(a,'Label',str);
    DATA.motorspeed = step./1000; %Custom is in uM/sec, not mm/sec
    DATA.customspeed = step./1000;
end
if DATA.motorspeed > 0
    ispeed = round(DATA.motorspeed.*DATA.speedscale.*DATA.stepscale/1000);
    fprintf('New Speed %.3f mm/sec  = %.0f RPM\n',DATA.motorspeed,ispeed);
    if ispeed < DATA.minrpm
        fprintf('%.1f less than min RPM (%.1f). Will use steps\n',ispeed, DATA.minrpm);
        fprintf(DATA.sport,'SP%.0f\n',DATA.minrpm);
    else
        fprintf(DATA.sport,'SP%.0f\n',ispeed);
    end
    DATA.motorspeed = ispeed;
end
SetMenuCheck(a,'exclusive');
set(DATA.toplevel,'UserData',DATA);


function v = rpm2speed(DATA)
v = DATA.motorspeed * 1000./(DATA.speedscale .* DATA.stepscale);

function CloseServoPort(a,b)
DATA = GetDataFromFig(a);
%tell verg first
if ~isempty(DATA.callback)
    feval(DATA.callback{:}, 'close', DATA);
end

if isfield(DATA, 'sport');
    fclose(DATA.sport);
    delete(DATA.sport);
end
close(DATA.toplevel);

function SetStepSize(a,b)

DATA = GetDataFromFig(a);
val = get(a,'value');
s = get(a,'string');
if strncmp(s(val,:),'Set Custom',8)
    step = DATA.customstep;
    defaultanswer{1} = num2str(step);
    str = inputdlg('Step Size (uM)','Custom',1,defaultanswer);
    if ~isempty(str)
    step = str2num(str{1});
    str = sprintf('%d (Custom)',step);
    DATA.customstep = step;
    DATA.step = step;
    s(val+1,1:length(str)) = str;
    set(a,'string',s,'value',val+1);
    end
else
    DATA.step = sscanf(s(val,:),'%d');
end
set(DATA.toplevel,'UserData',DATA);

function x = addfield(x,name, val)
if iscell(name) && iscell(val) && length(name) == length(val)
    for j = 1:length(name)
        x = addfield(x,name{j},val{j});
    end
elseif iscell(name)
    for j = 1:length(name)
        x = addfield(x,name{j},val);
    end
elseif ~isfield(x,name)
    x.(name) = val;
end

function DATA = SetDefaults(DATA)
DATA.position = 0;
DATA.step = 0;
DATA.stepscale = 65.6;
DATA.minrpm = 20;
DATA.speedscale = 15000;
DATA.timerange = 0;
DATA.firstsample = 0;
DATA.customspeed = 1; %also default
DATA.motorspeed = round(DATA.customspeed.* DATA.speedscale.*DATA.stepscale/1000);
DATA.setupfile = '/local/servomotor.setup';
DATA.logfile = '/local/servolog.mat';
DATA.lastspeed = 0;

DATA.motorid = -1; %< 0 means dont set id
DATA = addfield(DATA,{'stepsize' 'customstep' 'position'},0);
DATA = addfield(DATA,{'alldepths' 'alltimes' 'offidx'},[]);
txt = scanlines(DATA.setupfile,'silent');

if ~isempty(txt)
    if ~isfield(DATA,'ttyname')
        id = find(strncmp('serialport',txt,10));
        a = split(txt{id},'=');
        if isempty(id)
            DATA.ttyname = '/dev/tty.USA49Wfa1212P1.1';
        else
            DATA.ttyname = a{2};
        end
    end
    id = find(strncmp('position',txt,10));
    if ~isempty(id) && isempty(DATA.figpos) %if figpos set by command line, takes precedence
        DATA.figpos = sscanf(a{2},'%d');
    end
    id = find(strncmp('verbose',txt,10));
    if ~isempty(id)
        DATA.verbose = sscanf(a{2},'%d');
    end
elseif  ~isfield(DATA,'ttyname')
    cprintf('red','No Serial device Named, and missing file %s',DATA.setupfile);
end

if ~isfield(DATA,'figpos') || isempty(DATA.figpos)
    DATA.figpos = [1400 400 350 350];
end

for j = 1:length(txt)
    a = split(txt{j},'=');
    if strncmp(txt{j},'minrpm',6)
        DATA.minrpm = sscanf(a{2},'%f');
    end
end
if ~isfield(DATA,'callback')
    DATA.callback = [];
end
if ~isfield(DATA,'maxrange')
    DATA.maxrange = 3000000; %c. 30mm
end
if ~isfield(DATA,'plottype')
    DATA.plottype = 'FullHistory';
end

function SetHomePosition(DATA, pos)

if DATA.motorid >= 0
	fprintf(DATA.sport,sprintf('%dHO%.0f\n', DATA.motorid, pos .* DATA.stepscale));
else
	fprintf(DATA.sport,sprintf('HO%.0f\n', pos .* DATA.stepscale));
end
    set(DATA.setdepth,'string',sprintf('%.0f',pos));
    set(DATA.depthlabel,'string',sprintf('%.0f',pos));
    
    
function DATA = ChangePosition(DATA, step)

d = GetCurrentPosition(DATA,'show');
DATA = SetNewPosition(DATA,d+step);

function d = GetCurrentPosition(DATA, varargin)

showdepth = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'show',4)
        showdepth = 1;
    end
    j = j+1;
end

if DATA.motorid >= 0
fprintf(DATA.sport,sprintf('%dPOS\n',DATA.motorid));
else
fprintf(DATA.sport,sprintf('POS\n'));
end
%fprintf(DATA.sport,sprintf('0POS\n'));
pause(0.05);
[s, sdetails] = ReadLine(DATA.sport);
d = sscanf(s,'%d');
if DATA.verbose
    fprintf('From uDrive:%s\n',s);
end
if isempty(d)
    if ~isempty(s)
        fprintf('Microdrive response is not valid. Suggests configuation error in Serial Port.\n');
        fprintf('Microdrive Response was:%s\n',s);
        uiwait(warndlg(sprintf('MicroDrive Response Not Recognized. Serial Port Error?'),'Microdrive Error','modal'));
    elseif sdetails.nbtyes == 0
        uiwait(warndlg(sprintf('MircorDrive Serial Line No Response'),'Microdrive Error','modal'));                
    end
end
d = d./DATA.stepscale;
if showdepth
    set(DATA.depthlabel,'string',sprintf('%.0f uM',d));
end

function StopMotor(a,b)
DATA = GetDataFromFig(a);
if DATA.motorid >= 0
    fprintf(DATA.sport,'%dDI\n',DATA.motorid);
else
    fprintf(DATA.sport,'DI\n');
end
%If the speed was set to automatic, reset this
if DATA.motorspeed == 0 && DATA.lastspeed > 0
    DATA.motorspeed = DATA.lastspeed;
end


function DATA = SlowMove(DATA, pos)

newpos = pos .* DATA.stepscale;
d = GetCurrentPosition(DATA);
edur = abs(pos-d)./DATA.motorspeed;

nstep = floor(abs((pos-d)/10));  %divide into 10uM steps
smallstep = (pos -d)./nstep;
positions = d + smallstep .* [1:nstep];

if DATA.motorspeed > 0
    oldspeed = DATA.motorspeed;
elseif DATA.lastspeed > 0
    oldspeed = DATA.lastspeed;
else
    oldspeed = DATA.minrpm;
end

%real motor speed will be 100uM/sec, with 1 sec pauses
ispeed = round(0.1.*DATA.speedscale.*DATA.stepscale/1000);
fprintf(DATA.sport,'SP%.0f\n',ispeed);
DATA.motorspeed = ispeed;
stopb = findobj(DATA.toplevel,'tag','StopButton');
set(stopb,'enable','off');
stopui = findobj(DATA.toplevel,'tag','EmergencyStop');


DATA.firstsample = length(DATA.alldepths);
for j = 1:nstep
    istop = get(stopui,'value');
    if istop
        break;
    else
        DATA = SetNewPosition(DATA,positions(j));
        if j < nstep
            if DATA.motorspeed > 0 && edur < 1000
            pause(edur./nstep);
            else
                pause(1);
            end
        end
    end
end

if istop
    fprintf('Stopped by User\n');
    set(stopui,'value',0);
end
set(stopb,'enable','on');

fprintf(DATA.sport,'SP%.0f\n',oldspeed);
DATA.motorspeed = oldspeed;
DATA.firstsample = 0;

function DATA = SetNewPosition(DATA, pos)



newpos = pos .* DATA.stepscale;
pause(0.01);
d = NaN;

if DATA.motorspeed < DATA.minrpm;
    str = sprintf('%.1f less than min RPM (%.1f). Using steps',DATA.motorspeed, DATA.minrpm);
    fprintf('%s\n',str);
    set(gcf,'Name',str);
    DATA = SlowMove(DATA, pos);
    set(gcf,'Name','ServoMotor MicroDrive');
    return;
end

%first read position. 
if DATA.motorid >= 0
    fprintf(DATA.sport,'%dPOS\n',DATA.motorid);
else
    fprintf(DATA.sport,'POS\n');
end
    pause(0.01);
s = ReadLine(DATA.sport);
if isempty(s)
    uiwait(warndlg(sprintf('MicroDrive Not Responding'),'Microdrive Error','modal'));
    return;
end
if strcmp(s,'OK') %has not been sent ANSW1
    uiwait(warndlg(sprintf('MicroDrive Response not correct.\nTry "Reopen".'),'Microdrive Error','modal'));
    return;
end
d = sscanf(s,'%d');
if isempty(d)
    uiwait(warndlg(sprintf('MicroDrive Response Not Recognized'),'Microdrive Error','modal'));    
    return;
end
startpos = d./DATA.stepscale;

margin = 20 * DATA.stepscale; %allow for up to 20uM of noise
step = newpos-d;
edur = 0.1 + 2 .* abs(newpos-d) ./(DATA.motorspeed .* DATA.stepscale); %estimated duration
if edur > 600
    edur = 100;
end
if isempty(edur)
    edur = 1;
end
stopui = findobj(DATA.toplevel,'tag','EmergencyStop');
set(stopui,'value',0);
drawnow;

%set speed for movement
if DATA.motorspeed == 0 %automatic speed
    dd = abs(d-newpos)./DATA.stepscale; %step in uM
    if dd > 1000
        speed = 0.1;
    elseif dd > 50
        speed = 0.5;
    else
        speed = 1;
    end
    ispeed = round(speed.*DATA.stepscale.*DATA.speedscale/1000);
    if DATA.verbose
        fprintf('Speed %.3fmm/sec = %.0f rpm\n',speed,ispeed);
    end
    fprintf(DATA.sport,'SP%.0f\n',ispeed);
end

if DATA.verbose
    fprintf('Requesting %.0f->%.0f (~%.2f sec)\n',d,newpos,edur);
end

if abs(newpos) > DATA.maxrange
    uiwait(warndlg(sprintf('Requested Poition out of Range - Use Set button'),'Microdrive Error','modal'));
    return;
end


%request new position
if DATA.motorid >= 0
    fprintf(DATA.sport,'%dEN\n',DATA.motorid);
    fprintf(DATA.sport,sprintf('%dLA%.0f\n',DATA.motorid,newpos));
    pause(0.01);
    fprintf(DATA.sport,'%dM\n',DATA.motorid);
else
    fprintf(DATA.sport,'EN\n');
    fprintf(DATA.sport,sprintf('LA%.0f\n',newpos));
%    s = sprintf('LR%.0f\n',newpos-d)
%    fprintf(DATA.sport,s);
    pause(0.01);
    fprintf(DATA.sport,'M\n');
end
pause(0.01);
ts = now;
npost = 0;
newd(1) = d;

%Determine a range of acceptable postions given
%the step. So in case motor moves wrong way, or overshoots
%final position without ever having a small error, it is stopped.
startpos = d;
if step < 0
    minpos = startpos+step-margin;
    maxpos = startpos+margin;
else
    minpos = startpos-margin;
    maxpos = startpos+step+margin;
end


%keep reading postition during movement to check that all is well
%sometimes can stop just short of the finishing position, so allow some
%leeway.

j = 2;
while npost < 2
    if DATA.motorid >= 0
        fprintf(DATA.sport,sprintf('%dPOS\n',DATA.motorid));
    else
        fprintf(DATA.sport,sprintf('POS\n'));
    end
    s = ReadLine(DATA.sport);
    if strcmp(s,'OK') %in case Controller is in verbose state
        if DATA.verbose
            fprintf('Returned %s\n',s);
        end
        s = ReadLine(DATA.sport);
    end
    ts(j) = now;
    if ~isempty(s)
        try
            newd(j) = sscanf(s,'%d');
            set(DATA.depthlabel,'string',sprintf('%.0f uM',newd(j)./DATA.stepscale));
            drawnow;
            poserr = abs(newd(j)-newpos);
            if poserr < 3./DATA.stepscale.*1000 %3uM margin of error seems necessary. 
                npost = npost+1;
            end
            if edur > 2 || DATA.firstsample > 0
                PlotDepths(DATA, ts, newd);
            end
        catch
            newd(j) = 0;
            cprintf('error','Error reading Depth from %s: Verbose%d\n',s,DATA.verbose);
            fprintf('Returned %s\n',s);
        end
    else
        newd(j) = 0;
    end
    
%deatcitvat the motor if too much time, or close enough. 
    if mytoc(ts(1)) > edur+1
        fprintf(DATA.sport,'DI\n');
       npost = 2;
       F = gcf;
       if length(unique(newd)) == 1 %% no movment at all 
           uiwait(warndlg(sprintf('No movement reported at all.\n Check Cable to Microdrive',[newd(end) startpos newpos]./(DATA.stepscale.*1000)),'Movement incomplete Error','modal'));
       else
           uiwait(warndlg(sprintf('Failed to Complete Movement.\n Only Moved to %.3f (Requested %.3f->%.3f)',[newd(end) startpos newpos]./(DATA.stepscale.*1000)),'Movement incomplete Error','modal'));
       end
       figure(F);
    elseif newd(j) < minpos
        fprintf(DATA.sport,'DI\n');
       F = gcf;
       uiwait(warndlg(sprintf('Postion %.3f min allowd %.3f',[newd(end) minpos]./(DATA.stepscale.*1000)),'Microdrive request Error','modal'));
       figure(F);
    elseif newd(j) > maxpos
        fprintf(DATA.sport,'DI\n');
       F = gcf;
       uiwait(warndlg(sprintf('Postion %.3f Max allowed %.3',newd(end)./(DATA.stepscale.*1000),maxpos./(DATA.stepscale.*1000)),'Microdrive request  Error','modal'));
       figure(F);
    else
        stop = get(stopui,'value');
        if stop
            npost = 10;
        end
    end
    j = j+1;
end

DATA.position = newd(end)./DATA.stepscale;

dur = (ts(end)-ts(1)).*(24 * 60 *60);
dp = diff(newd);
id = find(abs(dp) > mean(abs(dp))/2); 
if length(id) > 1
    rate = (newd(id(end)+1) - newd(1))./(ts(id(end)+1) - ts(1));
    rate = rate ./ (DATA.stepscale * 1000 .* 60 * 60 * 24 );
else
    rate = (newd(end)-newd(1))/(dur.*DATA.stepscale.*1000);
end

if DATA.verbose
    fprintf('End pos %s took %.2f (%.2f) rate %.2f mm/sec\n',s,dur,edur,rate);
end

if DATA.motorid >= 0
    fprintf(DATA.sport,'%dDI\n',DATA.motorid);
else
    fprintf(DATA.sport,'DI\n');
end
pause(0.01);
DATA = PlotDepths(DATA, ts, newd);
SaveDiskLog(DATA);
it = findobj(DATA.toplevel,'tag','ManualSet');
if ~isempty(it)
    set(it,'string',sprintf('%.0f',DATA.position));
end

if ~isempty(DATA.callback)
    feval(DATA.callback{:}, newd(end)./DATA.stepscale, DATA);
end


function DATA = PlotDepths(DATA, ts, newd)

try
    
    if nargin == 1
        ts = DATA.newtimes;
        newd = DATA.newdepths;
    elseif nargin == 3
        DATA.newtimes = ts;
        DATA.newdepths = newd./DATA.stepscale;
        DATA.alltimes = [DATA.alltimes ts];
        DATA.alldepths = [DATA.alldepths DATA.newdepths];
        DATA.offidx(end+1) = length(DATA.alldepths);
    end
    DATA = CheckDepthData(DATA);
    if ~strcmp(DATA.plottype,'None')
        if strcmp(DATA.plottype,'LastMove')
            if DATA.firstsample > 0
                ts = DATA.alltimes(1+DATA.firstsample:end);
                y = DATA.alldepths(1+DATA.firstsample:end);
            else
                ts = ts-ts(1);
                y =newd./DATA.stepscale;
            end
            h = plot(ts, y);
            set(h,'ButtonDownFcn',@ServoPlotHit);
            
            set(gca,'xtick',[],'ytick',[],'ydir','reverse');
            xl = [min(ts) max(ts)];
            yl = [min(y) max(y)];
        elseif strcmp(DATA.plottype,'MoveSpeed')
            ts = ts-ts(1);
            k = 24 * 60 * 60 ./DATA.stepscale;
            dt = diff(ts) .* 24 * 60 * 60;
            y = diff(newd)./(dt .*DATA.stepscale);
            h = plot(ts(2:end), y);
            set(h,'ButtonDownFcn',@ServoPlotHit);
            
            set(gca,'xtick',[],'ytick',[],'ydir','normal');
            xl = [min(ts) max(ts)];
            yl = minmax(y);
            if length(y) > 3
                ms = prctile(y,50);
                line(xl,[ms ms],'color','r');
                text(mean(xl),ms,sprintf('%.1fuM/sec',ms),'verticalalignment','bottom');
            end
        else
            if DATA.timerange > 0
                ti = find(DATA.alltimes(end)-DATA.alltimes < DATA.timerange./(24 * 60));
                t = DATA.alltimes(ti);
                d = DATA.alldepths(ti);
            else
                t = DATA.alltimes;
                d = DATA.alldepths;
            end
            h = plot(t, d);
            set(h,'ButtonDownFcn',@ServoPlotHit);
            set(gca,'xtick',[],'ytick',[],'ydir','reverse');
            xl = minmax(t);
            yl = minmax(d);
        end
        
        if diff(yl) <= 0
            return;
        end
        tdur = diff(xl).*24; %hours
        tlabel = 'hr';
        if tdur < 1
            tdur = tdur .* 60;
            tlabel = 'min';
        end
        if tdur < 1
            tdur = tdur .* 60;
            tlabel = 'sec';
        end
        %N.B. ydir is reversed (deep == down)
        axis([xl yl]);
        text(xl(2),yl(1),sprintf('%.1f%s',tdur,tlabel),'horizontalalignment','right','verticalalignment','top');
        if strcmp(DATA.plottype,'MoveSpeed')
            text(xl(1),yl(2),sprintf('%.1fuM/sec',max(yl)),'horizontalalignment','left','verticalalignment','bottom');
            text(xl(1),yl(1),sprintf('%.1fuM/sec',min(yl)),'horizontalalignment','left','verticalalignment','top');
        else
            text(xl(1),yl(2),sprintf('%.1fuM',diff(yl)),'horizontalalignment','left','verticalalignment','bottom');
        end
        set(gca,'ButtonDownFcn',@ServoPlotHit);
    else
        delete(get(gca,'children'));
        bc = get(gcf,'color');
        set(gca,'color',bc,'box','off','xcolor',bc,'ycolor',bc);
    end
catch ME
    CheckExceptions(ME);
end

function ServoPlotHit(a,b)

pos = get(gca,'CurrentPoint');
cm = uicontextmenu;
uimenu(cm,'label',sprintf('%s %.1fuM',datestr(pos(1,1),'hh:mm'),pos(1,2)));
set(cm,'visible','on');



function [s, details] = ReadLine(port)

waitforbytes = 0.05; %wait 50ms if no bytes ready
ts = now;
while get(port,'BytesAvailable') == 0 && waitforbytes && mytoc(ts) < waitforbytes
    pause(0.001);
end
details.delay = mytoc(ts);
details.nbytes = get(port,'BytesAvailable');
s = fscanf(port,'%s');


function DATA = OpenServoPort(a,b)
 
DATA = GetDataFromFig(a);
x = instrfind('type','serial');
delete(x);
DATA.sport = serial(DATA.ttyname,'BaudRate',9600,'Timeout',1);
fopen(DATA.sport);
fprintf(DATA.sport,'ANSW1\n'); %?0 would stop unwanted "OK"
fprintf(DATA.sport,'SOR0\n'); %Source is Serial Line
fprintf(DATA.sport,'NET0\n');
fprintf(DATA.sport,'BAUD%d\n',9600);
fprintf(DATA.sport,'SP%d\n',DATA.motorspeed);

fprintf(DATA.sport,'APL1\n'); %ACtivate Position Limits
fprintf(DATA.sport,'LL%.0f\n',-DATA.maxrange);
fprintf(DATA.sport,'LL%.0f\n',DATA.maxrange);
fprintf(DATA.sport,'APL1\n');
if DATA.motorid >= 0
    fprintf(DATA.sport,'%dEN\n',DATA.motorid);
else
    fprintf(DATA.sport,'EN\n');
end


DATA.position = GetCurrentPosition(DATA,'show');
fprintf(DATA.sport,'DI\n');
set(DATA.toplevel,'UserData',DATA);
