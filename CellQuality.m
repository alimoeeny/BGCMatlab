function [Q, details] = CellQuality(C, varargin)
%[Q, details] = CellQuality(C, varargin) generates a single number for 
% goodness of a cell, including isolation and drop stats
% Q  = 1./(1./(dropi-2) + 0.75./(isol-1));
% so, Q > 1 
% CellQuality(Dropi, Isol)  calculates value for a set of isolaiton/drop
% intex numbers
isolation = 'bestisolation';
details.isolation_mode = isolation;
c = [];

if strcmp(C,'selfplot')
    TestPlot(varargin{:});
    return;
end
j = 1;
while j <= length(varargin)
    if isnumeric(varargin{j})
        c = varargin{j};
    end
    j = j+1;
end

%if called with CellQuality(Matrix1, Matrix2)
if isnumeric(C) && sum(size(c) == size(C)) == 2
    d = C;
end
if isempty(c)
    if strcmp(isolation,'bestisolation')
        if isfield(C,'bestisolation')
            if isfield(C.bestisolation,'isolation')
                d = C.bestisolation.isolation(1);
            elseif isfield(C,'isolation');
                d = C.isolation(1);
            end
        elseif isfield(C,'myisolation')
            d = C.myisolation(1);
        elseif isfield(C,'isolation')
            d = C.isolation(1);
        else
            d = NaN;
        end
    end

    if isfield(C,'dropi') && length(C.dropi) > 2
        c = C.dropi(3);
    else
        c = NaN;
    end
end
%dropi is -1 ->1, isolation is -3 ->inf, but -3 not really poss

dth = 2; %for isoltion
cth = 1;
c(c<cth) = cth;
d(d<dth) = dth;
Q = erf((c-2)./2) + (d-3);
Q = 0.75./(c-cth) + 1./(d-dth);
Q = 1./Q;


function TestPlot(varargin)

dx = [0:0.1:10];
[D,I]  = meshgrid([0:0.1:10],[0:0.1:10]);
im = CellQuality(D,I);
subplot(2,1,1);
imagesc(dx,dx,im);
colorbar;
subplot(2,1,2);
hold off;
plot(dx,im(end,:));
hold on;
plot(dx,im(:,end));
legend('D','I');

