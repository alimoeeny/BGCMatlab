function varargout =  mymemory(varargin)
%[used, physicalsize] = mymemory(...) wrapper for memory that avoids crashes
%

try 
    [x,y] = memory();
    varargout{1} = x.MemUsedMATLAB;
    varargout{2} = y.PhysicalMemory.Total;
    varargout{3} = y.PhysicalMemory.Available;
catch
    varargout{1} = NaN;
    varargout{2} = NaN;
end