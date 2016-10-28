    function name = NumberToExptName(x, Expt)
        name = 'e0';
        if isempty(x)
            return;
        end
        switch x
            case 107
            case 112.5
            name = 'sf';
            case 101
            name = 'me';
            case 96
            name = 'Pp';
            case 106
            case 100
                case 250
            name = 'Op';
            case 111
                name = 'or';
            case 6
                name = 'dx';
            case 117
                name = 'ce';
            case 9
                name = 'dx';
            case 119
            case 237
            case 132
                name = 'ce';
            case 6
                name = 'dx';
            otherwise
            fprintf('Unknown');
        end
        
        fprintf('Ex %d ->%s\n',x,name);


%    etvals = [6 117  9 132 119 237 250 111 250 7];
%    fprintf('%s Expts %dX%d\n',name,Expt.Stimvals.et,Expt.Stimvals.e2);
%    etnames = {'dx' 'ce' 'dx' 'ce' 'ce' 'dp' 'Op' 'e0'};