function S = CellToStruct(C, varargin)
% S = CellToStruct(C, varargin)
%converts a cell array to a structure. N.B quite different from cell2struct
%This is Inteded for cell arrays that have similar elements. cell2struct
%takes all the fields in each element of C, and then makes a structure with
%these fields.
%Empty fields are set to NaN, so that cat(S.field) yeilds N elsemnt
%S = CellToStruct(C, 'nopad') leaves empty fields empty
nanpad = 1;
flatcat = 0;
strarg = cell2cellstr(varargin);
if sum(strcmp('nopad',strarg))
    nanpad = 0;
end
if sum(strcmp('flat',strarg))
    flatcat = 1;
end
sz = size(C);
C = C(:);
S = [];

allfields = {};
nr = 0;
for j = 1:length(C);
    if isstruct(C{j})
        nr = nr+1;
        f = fields(C{j});
        allfields = unique({f{:} allfields{:}});
        for s = 1:length(C{j})
            if flatcat
                ns = length(S)+1;
                nt  = 1;
            else
                ns = s;
                nt = nr;
            end
            for k = 1:length(f)
                S(nt,ns).(f{k}) = C{j}(s).(f{k});
            end
        end
    end
end
if ~exist('S','var')
    S = [];
    return;
end

if nanpad
    for j = 1:length(S)
        for k = 1:length(allfields)
            f = allfields{k};
            if isempty(S(j).(f))
                if iscell(S(j).(f))
                    S(j).(f) = {};
                else
                    S(j).(f) = NaN;
                end
            end
        end
    end
end
if size(S,2) == 1 && flatcat == 0 && sz(2) > 1 %make rows/vols match
    if prod(size(S)) == prod(sz) %no mixup caused by empty cells
       S = reshape(S,sz);
    else
       sz(2) = round(prod(size(S))./sz(1)); 
       S = reshape(S,sz);
    end
end