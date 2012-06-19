function [hostname, details] = gethostname(varargin)

os = computer;
if strncmp(os,'PCWIN',5)
    hostname=getenv('USERDOMAIN');
end
details.hostid = hostid;
