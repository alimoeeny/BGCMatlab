function [X, xid, details] = FindXCorr(xcorrs, id)
%[X, xid] = PC.FindXCorr(xcorrs, id) find element of xcorrs struct
% with data for probe pair id.  
%id can either be floats defining P/Cl or ids to xcorrs.p
%xcorrs can also be a struct containing xcorrs (e.g. CellGroup)
%[X, id, details] = PC.FindXCorr(xcorrs, id)
%               details.efficacy has efficacy 5/6 (adjusted scores
%               but ordered so that efficacy 1 is alwasy id(1)->id(2)
if isfield(xcorrs,'xcorrs')
    [X,xid, details] = PC.FindXCorr(xcorrs.xcorrs,id);
    if isfield(xcorrs,'eid')
        X.eid = deal(xcorrs.eid);
    end
    return;
end
if sum(rem(id,1)) > 0
    p = cat(1,xcorrs.probe);
else
    p = cat(1,xcorrs.p);
end
xid = find((p(:,1) == id(1) & p(:,2) == id(2)) | (p(:,1) == id(2) & p(:,2) == id(1)));
X = xcorrs(xid);
if nargout > 2
    if p(xid,1) == id(1)
        details.efficacy = X.efficacy(5:6);
        details.nspk = X.nspk;
    else
        details.efficacy = X.efficacy([6 5]);
        details.nspk = X.nspk([2 1]);
    end        
end