function FitDriftMeanSpike(M, varargin)

npts = 40;
mx = 8.*npts;
T = MakeFirstMean(M,0, mx);
for j = 1:10;
    pos = GetAlignment(M,T, mx);
    T = MakeAlignedMean(M,pos,mx);
    offset = minmax(pos);
    allpos(j,:)  = pos;
    AllT{j} = T;
end
imagesc(allpos);


function T = MakeFirstMean(M,pos,mx)

nprobes = 24;
npts = 40;
T = ones(size(M,1),size(M,2)+diff(minmax(pos)).*npts).*NaN;
T = M(end,:);
l = size(M,2);
AllT(1,mx:mx+length(T)-1) = T;
AllT(1,1:mx-1) = NaN;
AllT(1,end:end+mx) = NaN;
for j = size(M,1)-1:-1:1
    x(j,:) = ShiftXcorr(nanmean(AllT,1), M(j,:),npts, nprobes,mx);
    [a,b] = max(x(j,:));
    b = b-8;
    if b.*npts > -mx
        AllT(j,:) = NaN;
        o = mx+b.*npts;
        AllT(j,o:o+l-1) = M(j,:);
    end
    pos(j) = b;
end
T = nanmean(AllT,1);


function T = MakeAlignedMean(M,pos,mx)

npts = 40;
T = ones(size(M,1),size(M,2)+2.*mx).*NaN;
for j = 1:size(M,1)
    apts = [1:size(M,2)] + (pos(j).*npts) + mx;
    T(j,apts) = M(j,:);
end
T = nanmean(T);

function [pos, offset] = GetAlignment(M,T, offset)

fnpts = 40;
nprobes = 24;
for j = 1:size(M,1);
        x(j,:) = ShiftXcorr(T, M(j,:),fnpts, nprobes,offset);
        [a, b] =max(x(j,:));
        mxcs(j) = b-8;
        mxval(j) = a;
end
offset = [0 0];
pos = mxcs;

function xcs = ShiftXcorr(allshape, rowshape, npts, nprobes, offset)
    xpts = (length(allshape) - npts .*nprobes);
    nprobes = round(length(rowshape)./npts);
    for k = [-7:7]+offset
        if k > 0
            bpts = 1:length(rowshape);
            apts = k:k+length(bpts)-1;
        end
        id = find(~isnan(allshape(apts)));
        apts = apts(id);
        bpts = bpts(id);
        xc = corrcoef(allshape(apts),rowshape(bpts));
        xcs(k+8-offset) = xc(1,2);
    end
