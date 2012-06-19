os = computer;
if strmatch(os,{'MAC' 'MACI' 'GLNX64' 'MACI64'})
path('/bgc/bgc/matlab/dev',path);
elseif exist('Z:/bgc/matlab','dir')
path('Z:/bgc/matlab/dev',path);
else
path('C:/bgc/bgc/matlab/dev',path);
end