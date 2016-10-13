function took = mytoc(tstart, tend)
%took mytoc(start)
%calculates elapsed time since start (start = now)
%if no return value is requested, prints result.
if nargin == 1
    tend = now;
end
took = (tend-tstart) * 24 * 60 * 60;
if nargout == 0
   fprintf('Took %.3f sec\n',took);
end
       
