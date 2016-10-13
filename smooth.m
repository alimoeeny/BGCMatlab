function smoothed = smooth(x, w, varargin)
%smooth(x, w, ...) smooth data with boxcar, width w
%  if w is omitted, 5 is used. 
%smooth(x,w,'gauss') uses Gaussian smoothing with an SD of w
%if x is a matrrix, just smoohthsalong dim 2;
% smooth(x, w, 'dimension',2) smooths along 1. 

if nargin < 2
    w = 5
end

sdim = 1;
kernel = [];
j = 1;
if strncmpi(w,'kernel',4)
    w = 1;
    kernel = varargin{1};
    j = 2;
end

if iscell(x)
    smoothed = x; %to set dimenstions
    for j = 1:length(x(:))
        smoothed{j} = smooth(x{j},w, varargin{:});
        smoothed{j} = reshape(smoothed{j},size(x{j}));
    end
    return;
end
while j < nargin-1
    if strncmpi(varargin{j},'gauss',4)
        if w > floor(length(x)./4);
            w = floor(length(x)./4);
        end
        kernel = Gauss(w,[-(w*2):(w*2)]);
        kernel = kernel./sum(kernel);
    elseif strncmpi(varargin{j},'dimension',3)
        j = j+1;
        sdim = varargin{j};
    elseif strncmpi(varargin{j},'kernel',4)
        j = j+1;
        kernel = varargin{j};
    end
    j = j + 1;
end
if w == 0
    smoothed = x;
    return;
end
if isempty(kernel)
    kernel = ones(1,w)/w;
end
if length(x) < length(kernel)
    smoothed = ones(size(x)).* mean(x);
    return;
end
if numel(x) >  length(x) % not a vector
    if sdim == 1
        smoothed = conv2(kernel,sdim,squeeze(x),'same');
    else
        smoothed = conv2(sdim,kernel,squeeze(x),'same');
    end
else
smoothed = conv(x,kernel);
start = ceil(length(kernel)/2);
last = start + length(x) - 1;
smoothed = smoothed(start:last);
for j = 1:start
    smoothed(j) = smoothed(j) .* sum(kernel)/sum(kernel(1:(start+j-1)));
    smoothed(end-j+1) = smoothed(end-j+1) .* sum(kernel)/sum(kernel((start-j+1):end));
end
end