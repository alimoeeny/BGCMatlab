function d = int2double(x, scale, varargin)
%int2double(x, scale) convert in to double, mapping 32002 ->NaN

maxint = double(intmax('int16')-5);
d = double(x).* scale(1) ./maxint;
d(x==maxint+2) = NaN;
