function pos = MakeTriplet(X, varargin)
%PC.MakeTripet converts strutures into 3xN matrices [e p c; ...]
%for easy access to cluster info from EG. a CellGroup
pos= [];

if isfield(X,'groups') && isfield(X,'cells') %A Duplcicates listing for an Expt
    if isnumeric(varargin{1})
        id = varargin{1};
    else
        id = 1; %Dummy for now
    end
    e = X.eid;
    for j = 1:length(id)
        p = floor(X.xcorrs(id(j)).probe);
        cl = round(rem(X.xcorrs(id(j)).probe,1).*10);
        pos(end+1,:) = [e p(1) cl(1)];
        pos(end+1,:) = [e p(2) cl(2)];
    end
end