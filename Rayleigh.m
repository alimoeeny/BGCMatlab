function [r, pval, theta] = Rayleigh(angles, varargin)
% [r, pval, theta] = Rayleigh(angles, varargin)
% comput vertor mean of a set of angles
nresample = 0;
axialdata = 0;
pval = NaN;
j = 1;
while j < nargin
    if strncmpi('axial',varargin{j},3)
        axialdata =1;
    elseif strncmpi('resample',varargin{j},3)
        j = j+1;
        nresample = varargin{j};
    end
    j = j+1;
end

if axialdata
    angles = angles .* 2;
end

sina = mean(sin(angles));
cosa = mean(cos(angles));
r = abs(cosa + i * sina);
theta = atan2(sina,cosa);

if nresample
   rnd = rand(length(angles),nresample) .* 2 * pi;
   rs = abs(i * mean(sin(rnd)) + mean(cos(rnd)));
   pval = sum(rs >= r)./nresample;
end

if axialdata
    theta = theta./2;
end
