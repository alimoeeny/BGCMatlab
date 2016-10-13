ver
ts = now;
os = computer;
if exist('homepath','var') %path already set
elseif strmatch(os,{'MAC' 'MACI' 'MACI64'})
    path(path,'/b/group/matlab/agb');
    path(path,'/b/group/matlab');
    path('/b/bgc/matlab',path);
%    path('/b/bgc/matlab/dev',path);
    cd /b/bgc/matlab
elseif strmatch(os,'GLNXA64')
    path(path,'/b/group/matlab');
    path('/b/bgc/matlab',path);
    cd /b/bgc/anal;
elseif strmatch(os,{'PCWIN' 'PCWIN64'})
    if exist('Z:/bgc/matlab','dir')
        path('Z:/bgc/matlab',path);
        path(path,'Z:/group/matlab');
        cd Z:/bgc/matlab;
    else
        path('/b/bgc/matlab',path);
        path(path,'/b/group/matlab');
        cd /b/bgc/matlab;
    end
end

hostname = gethostname;
hostname = strrep(hostname,'.nei.nih.gov','');
hostname = strrep(hostname,'.','_');
hostfile = ['startup_' hostname];
if exist(hostfile)
    fprintf('Running %s for host specific settings\n',hostfile);
    eval(hostfile);
end
precompile = 0; %set this when testing startup
if strcmp(hostname,'NEIK2C36BC13P') && precompile
    fprintf('initiating src cache for sucombine\n');
    sucmb.InitSrcCache;
end
fprintf('Starting with %s under %s took %.2f on %s\n',mfilename('fullpath'),os,mytoc(ts),datestr(now));
dbstop if error;

