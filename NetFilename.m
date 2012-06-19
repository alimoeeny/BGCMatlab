function [netname, netdir] = NetFilename(name, varargin)

[a,b] = fileparts(name);

if strfind(name,'/lem')
    monkey = 'lem';
elseif strfind(name,'/jbe')
    monkey = 'jbe';
else
    monkey = 'lem';
end

netname = regexprep(name, ['.*/' monkey '/'],['Z:/bgc/data/' monkey '/']);
if strcmp(netname, name) %no match
    netname = [];
    netdir = 0;
else
    [a,b] = fileparts(netname);
    netdir = a;
end