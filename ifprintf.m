function ifprintf(DATA, varargin)
%ifprintf(DATA, varargin{:}) call fprintf if DATA.profiling


str = sprintf(varargin{:});
if isfield(DATA,'firstcalltime') && DATA.firstcalltime > 0
    ts = sprintf('at %.3f', mytoc(DATA.firstcalltime));
    if regexp(str,'\n$');
        str = regexprep(str,'\n$',[ts '\n'])
    else
        str = [str ts];
    end
end
if DATA.profiling
    fprintf('%s',str);
end