function distance = DistanceToEllipse(E, pos, varargin)   strargs = cell2cellstr(varargin);%When called from a new graph, ther is nothing to check, but E may have%leftover bits inif isempty(E) | ~isfield(E,'pos') | ~isfield(E,'angle');    if isfield(E,'ellipse') && isfield(E.ellipse,'xy')        X.xyr = E.ellipse.xy(1:4);        X.angle = E.ellipse.xy(5);        X.shape = 0;        X.pos = X.xyr(1:2);        distance = AllV.DistanceToEllipse(X, pos, varargin(:));        return;    end    distance = [NaN NaN];    return;endif E.shape == 1    r(1) = AllV.LineLength(E.pos)/3;    r(2) = r(1)./10;    a(1) = (E.pos(3)+E.pos(1))/2; %x radius    a(2) = (E.pos(4)+E.pos(2))/2;    xy = pos - a;    xy = xy ./r;    cn = cos(-E.angle);    sn = sin(-E.angle);    p(1) = xy(1) * cn + xy(2) * sn;    p(2) = xy(2) * cn - xy(1) * sn;    distance = sum(p.^2);    distance(2) = NaN;else            if isfield(E,'xyr') && length(E.xyr) > 3                if length(E.xyr) > 3                    r = E.xyr(3:4);                else                    r = 0;                end                a = E.xyr(1:2);            else                r(1) = (E.pos(3)-E.pos(1))/2; %x radius                r(2) = (E.pos(4)-E.pos(2))/2;                a(1) = (E.pos(3)+E.pos(1))/2; %x radius                a(2) = (E.pos(4)+E.pos(2))/2;            end        xy = pos - a;        xy = xy ./r;        cn = cos(-E.angle);        sn = sin(-E.angle);        p(1) = xy(1) * cn + xy(2) * sn;        p(2) = xy(2) * cn - xy(1) * sn;%        p = p./r;        distance(1) = sum(p.^2);    if isfield(E,'markpt')        distance(2) = sqrt(sum((pos -E.markpt).^2))./mean(r);    else        distance(2) = 1;    endendif sum(strcmp('radius',strargs))                distance = distance(1);end