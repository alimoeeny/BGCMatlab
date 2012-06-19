function [monk, monkyename] = GetMonkeyName(name)
%[monk, monkyename] = GetMonkeyName(name)
%

if findstr(name,'lem')
    monk = 'lem';
    monkeyname = 'lem';
elseif findstr(name,'dae')
    monk = 'dae';
    monkeyname = 'dae';
elseif findstr(name,'jbe')
    monk = 'jbe';
    monkeyname = 'jbe';
else
    fprintf('Can''t find monkey name is %s\n',name)
    monk = 'nby';
    monkeyname = 'none';    
end