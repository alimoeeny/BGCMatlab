function FullV = fix(FullV, varargin)
%fullv.fix(FullV, 'toint') converts double to int
%fullv.fix(FullV, 'todouble') converts int to double

strargs = {};
j= 1;
while j <= length(varargin)
    if ischar(varargin{j})
        strargs = {strargs{:} varargin{j}};
    end
    j = j+1;
end


if sum(strcmp('todouble',strargs)) && isinteger(FullV.V)
    if isfield(FullV,'intscale')
    FullV.V = double(FullV.V);
    np = size(FullV.V,1);
    for j = 1:np
        FullV.V(j,:) = FullV.V(j,:) .* FullV.intscale(1)/FullV.intscale(2);
    end
    end
end
