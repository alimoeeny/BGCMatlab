function [d, scale] = double2int(x, scale, varargin)
%[i, scale] = double2int(x, scale) convert double to int , mapping NaN ->32764 
%[i, scale] = double2int(x) determine scale from max(abs(x(:)))

if nargin == 1
    scale = max(abs(x(:)));
end

maxint = double(intmax('int16')-5);
d = int16(round(x.*maxint./scale));;
d(isnan(x)) = maxint+2;
