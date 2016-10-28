function Expt = CRData2int(Expt, varargin)
%Called by spike2 to store data as int instead of double

if isfield(Expt,'Trials') && isfield(Expt.Trials,'Eyevals')
        for j = 1:length(Expt.Trials)
            f = fields(Expt.Trials(j).Eyevals);
            for k = 1:length(f)
                x = minmax(Expt.Trials(j).Eyevals.(f{k}));
                if length(x) == 2
                    allx(k,:) = x;
                end
            end
            x = minmax(allx);
            if length(x) == 2
                emrange(j,:) = x;
            end
        end
        emrange = minmax(emrange(:));
        emscale =  max(abs(emrange));
        Expt.Header.emscale(1) = emscale;
        maxint = double(intmax('int16')-5);
        Expt.Header.emscale(2) = maxint;        
        else
            emscale = Expt.Header.emscale(1);
            maxint = Expt.Header.emscale(2);
        end
        for j = 1:length(Expt.Trials)
            f = fields(Expt.Trials(j).Eyevals);
            for k = 1:length(f)
                x = Expt.Trials(j).Eyevals.(f{k});
                x = round(x .*maxint./emscale);
                x(isnan(x)) = maxint+2;
                Expt.Trials(j).Eyevals.(f{k}) = int16(x);
            end
        end
        Expt.Header.EMdataclass = 'int';    
end