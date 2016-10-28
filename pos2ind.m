function [e,p,cl] = pos2ind(pos, varargin)
%[e,p,cl] = pos2ind(pos, varargin) convert float E-P into E,P,cl
%for clusters etc
%if pos is scalar returns probe,cluster,c (also = expt, subexpt,c);
%       where c = 'a' if subexpt = 1, 'b' for 1 etc
%else returns expt,probe,cluster

if length(pos) == 1
    p = round(rem(pos,1).*10);
    e = floor(pos);
    cl = char('a' +p);
else
    cl = round(rem(pos(2),1).*10);
    p = floor(pos(2));
    e = floor(pos(1));
end