function C = AddFits(DATA, C, varargin)
%PC.AddFits(DATA,C)
%Adds back in fit data to the Clusters in cell array C

if isfigure(DATA)
    F = DATA;
else
    F = DATA.toplevel;
end

CF = getappdata(F,'AutoClusterFits');
if iscell(C)
    for j = 1:length(C)
        if iscell(C{j})
            C{j} = PC.AddFits(DATA, C{j}, varargin{:});
        else
        if isfield(C{j},'probe')
            p = C{j}.probe;
        else
            p = j;
        end
        X = CF{C{j}.exptid}{p};
        C{j} = CopyFields(C{j},X);
        end
    end
else
    X = CF{C.exptid}{C.probe};
    C = CopyFields(C,X);
end
   

