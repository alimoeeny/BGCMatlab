function latency = sdflatency(sdf, presamples,offset)
%
% sdflatency(sdf, presamples)
% returns a latency estimate from a spike density function, sdf
% presamples is the number of samples in the sdf guaranteed not to
% contain a response.
% offset, if defined, determines a starting location in the buffer

if nargin < 3
  offset = 1;
end

if presamples+offset > length(sdf)
    latency = NaN;
    return;
end
prerate = mean(sdf(offset:presamples+offset));
ti = presamples+offset;
count = 0;
while ti < length(sdf) & count < 5
  if sdf(ti) > prerate
    count = count + 1;
  else
    count = 0;
  end
  ti = ti + 1;
end
latency = ti-count;
