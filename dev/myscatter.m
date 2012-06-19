function myscatter(x,y,symbol,varargin)
%myscatter(x,y,symbol, ... makes a scatterplot where each datapoint has a callback
%function. Default is to print out index of data point in the vectors x,y
%
% myscatter(x,y,'.','buttondown',@fcn)
%               sets the function called after a press
% myscatter(x,y,'.','ids',idlist)
%               sets the id numbers that are passed to @fcn


fcn = @scatterhit;
idlist = [];
plotargs = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'buttonpress',10)
        j = j+1;
        fcn = varargin{j};
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        idlist = varargin{j};
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end

if isempty(idlist)
    idlist = 1:size(x,1);
    bidlist = 1:size(x,2);
end

for k = 1:size(x,1)
for j = 1:size(x,2)
    plot(x(k,j),y(k,j),symbol,'buttondownfcn',{fcn, idlist(k), bidlist(j)},plotargs{:});
    hold on;
end
end


function scatterhit(a,b,ida, idb)

fprintf('Point %d, %d\n',ida, idb);

