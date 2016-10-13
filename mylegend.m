function lh = mylegend(h, labels, varargin)
% mylegend(h, labels, varargin)
% calls legend with a list of handles and labels, but checks these
% first to avoid annoying legend crashes.
ilb = []; %in case labels is empty.
if isempty(labels)
    return;
end
lh = [];
%in case theses are matrices, collapse to 1D.
legendloc = '';
skip = [];
j =1;
while j <=length(varargin)
    if sum(strcmpi(varargin{j},{'LowerRight' 'BottomOutside'}))
        legendloc = varargin{j};
        skip(end+1) = j;
    end
    j = j+1;
end
varargin = varargin(setdiff(1:j-1,skip));



h = h(:);
labels = labels(:);
for j = 1:length(labels)
    if ischar(labels{j})
        ilb(j) = 1;
    else
        ilb(j) = 0;
    end
end
id = intersect(find(myhandle(h)),find(ilb));
if length(id)
    try
        lh = legend(h(id),{labels{id}}, varargin{:});
        if ~isempty(legendloc)
            gui.PlaceLegend(lh, legendloc);
        end
    catch ME
        CheckExceptions(ME);
    end
end

