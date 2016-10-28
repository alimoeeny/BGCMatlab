function F = KeepFigure(a, varargin)
%KeepFigure(figure,  ...) Changes a figure Tag and name  (prepending 'Keep') 
%so that subsequent calls to GetFigure do NOT
%return this figure.  
%figure can be a figure handle or a Tag.  
%KeepFigure(figure,  'detach'...) reomves any association with a parent, so
%KeepFigure(figure,  'label', string) prepends with string instead of 'Keep'
%that Closing the the parent will not close the kept figure

lbl = 'Keep';
if nargin < 1
    a = gcf;
end
if isfigure(a)
    F = a;
else
    F= GetFigure(a);
end

replaceduplicates = 0;
detatch = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'detatch',4)
        detatch = 1;
    elseif strncmpi(varargin{j},'label',4)
        j = j+1;
        lbl = varargin{j};
    elseif strncmpi(varargin{j},'replace',4)
        replaceduplicates = 1;
    end
    j = j+1;
end

tag = get(F,'Tag');
newtag = [lbl tag];
dpl = findobj(allchild(0),'flat','type','figure','tag',newtag);
if  ~isempty(dpl)
    if replaceduplicates
        fprintf('KeepFigure:  Deleting existing figure tag %s\n',newtag);
        close(dpl);
    else
        fprintf('KeepFigure:  Duplicates figure tag %s\n',newtag);
    end
end
 
set(F,'tag', [lbl tag]);
name = get(F,'name');
set(F,'name', [lbl name]);
setappdata(F,'KeepFigure',tag); %record original tag
if detatch
    if isappdata(F,'ParentFigure')
        rmappdata(F,'ParentFigure');
    end
end