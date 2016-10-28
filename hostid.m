function r = hostid()
%my  function to return host name
[~, r] = system('hostname');