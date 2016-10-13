function CloseTag(tag)
it = findobj(allchild(0),'flat','Tag',tag);
for j = 1:length(it)
  fcn =   get(it(j),'CloseRequestFcn');
  if isempty(fcn) || ischar(fcn)
      delete(it(j));
  else
      fcn{1}(it(j), fcn{:});
  end
end

