function [hostname, details] = gethostname(varargin)

os = computer;
if strncmp(os,'PCWIN',5)
    hostname=getenv('USERDOMAIN');
else
    [a, hostname] =system('hostname');
end
details.hostid = hostid;
