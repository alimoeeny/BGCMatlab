function res = parexec(cmds, varargin)
%res = parexec(cmds) calls eval with each string in the cellstr cmds, in a parfor
%loop
res = [];

for j = (1:length(cmds))
    try
        eval(cmds{j});
    end
end