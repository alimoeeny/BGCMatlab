function dur = save(FullV, varargin)
%fullv.save(FullV) saves fullV to FullV.loadname
%converts double back to int first
%adds -v7.3 if size > 2Gb.
verbose = 0;
strargs = cell2cellstr(varargin);
if sum(strcmp('verbose',strargs))
    verbose = 1;
    id = find(strcmp('verbose',strargs));
    varargin = varargin(setdiff(1:length(varargin),id));
end

if isfloat(FullV.V)
    intscale(2) = double(intmax('int16')-5);
    intscale(1) = max(abs(FullV.V(:)));
    FullV.V = int16(round(double(FullV.V .* intscale(2)./intscale(1))));
    FullV.intscale = intscale;
end
x = whos('FullV');
ts = now;
if verbose
    fprintf('Saving %.2f Gb to %s...',x.bytes./(1024.*1024 .* 1024),FullV.loadname);
end
if x.bytes > 2.1e+09
    save(FullV.loadname,'FullV',varargin{:},'-v7.3');
else
    save(FullV.loadname,'FullV',varargin{:});
end
dur=mytoc(ts);
if verbose
    fprintf('took %.2f\n',dur);
end
