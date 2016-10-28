function [estb0,estb1,varargout] = fit_bothsubj2error(x,y,varargin)
% function [estb0,estb1,varargout] = fit_bothsubj2error(x,y,varargin)
% This is a function implementing Draper & Smith's, p.91, recipe for fitting a straight line between variables x and y,
% where *both* variables are subject to error.
% lambda is the ratio of the variance on y to the variance on x. If unspecified, it is taken to be 1.
% In the limit lambda=infinity, there is no error on the x-values, and this algorithm becomes the standard linear regression.
% I find its calling syntax easier than that of Matlab's regress function.
% The fitted line is y = estb0 + estb1 x
% Optional output arguments can give you the 95% confidence interval on estb0 and estb1,
% which are obtained by resampling. An optional 4th argument gives the number of
% resampling runs to use (the default is 1000)
%
% Example: 
% x=rand(1,100); y= rand(1,100); lambda = 2;
% [estb0,estb1] = fit_bothsubj2error(x,y)
% [estb0,estb1] = fit_bothsubj2error(x,y,lambda)
% [estb0,estb1,b0CI,b1CI] = fit_bothsubj2error(x,y)
% [estb0,estb1,b0CI,b1CI] = fit_bothsubj2error(x,y,lambda,5000)
% NB:
% [estb0,estb1] = fit_bothsubj2error(x,y,Inf) gives the same answer as Matlab's own b = regress(y',[ones(size(x));x]')
% with a more convenient syntax.
%
% Jenny Read 2/4/2003

if length(x)~=length(y)
    disp('I need x and y to be the same length!')
    return
end

if nargin>2
    lambda = varargin{1};
else
    lambda = 1;
end

if lambda<1
    disp('Warning!! Do not try and resample with lambda<1. It is a better idea to swap your data over and get lambda > 1')
end


yWeightArg = strncmp(varargin,'yweights',4);
if any(yWeightArg)
    yWeights = varargin{find(yWeightArg)+1};
    ybar = wmean(y,yWeights);
    xbar = wmean(x,yWeights);
else
    ybar = mean(y);
    xbar = mean(x);
    yWeights = ones(size(y));
end
Syy = sum((y-ybar).^2);
Sxx = sum((x-xbar).^2 .* yWeights);
Sxy = sum((x-xbar).*(y-ybar).*(yWeights));

if isinf(lambda)
    % then the variance on x is negligible in comparison to the variaance on y; go back to the standard regression
    estb1 = Sxy / Sxx;
elseif lambda==0
    estb1 = Syy / Sxy;
else
    estb1 = ( Syy - lambda * Sxx + sqrt((Syy-lambda*Sxx).^2 + 4*lambda*Sxy^2) ) / (2*Sxy);
end
estb0 = ybar - estb1 * xbar;

if nargout==4
  % Resample
  if nargin>3
      nresampleruns = varargin{2};
  else
      nresampleruns = 1000;
  end
   n = length(x);
   for jrun=1:nresampleruns
       indx = unidrnd(n,1,n);
       newx = x(indx);
       newy = y(indx);
       % call the routine recursively!
       [newestb0(jrun),newestb1(jrun)] = fit_bothsubj2error(newx,newy,lambda,nresampleruns);
   end
   % Output confidence intervals for b0, b1
   varargout{1} = [prctile(newestb0,2.5) prctile(newestb0,97.5)];
   varargout{2} = [prctile(newestb1,2.5) prctile(newestb1,97.5)];
   
end