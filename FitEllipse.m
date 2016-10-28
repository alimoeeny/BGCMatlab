function ell = FitEllipse(xy, varargin)
%given a set of xy pairs. Find the neareast ellipse that matches
%the shape, by minimizing distance from pts to teh ellipse boundary
guess = [];
aspectratio = 1;
j=1;
while j <= length(varargin)
    if strcmp(varargin{j},'aspectratio')
        j = j+1;
        aspectratio = varargin{j};
    elseif strcmp(varargin{j},'show')
        ShowFits(xy, aspectratio);
        return;
    end
    j= j+1;
end
if isempty(guess)
    angles = 0:0.1:pi;
    for j = 1:length(angles)
        xyr([1 2]) = mean(xy);
        xyr(5) = angles(j); %angl
        xys = xyrotate(xy(:,1)-xyr(1),xy(:,2)-xyr(2),xyr(5));
        xyr([3 4]) = std(xys).*1.4;
        scores(j) = TryEllipse(xyr,  xy, aspectratio);
        xyrs(j,:) = xyr;
    end
    [a,b] = min(scores);
    xyr = xyrs(b,:);
else
end
trackfit = 1;
if trackfit
    [score,rs] = TryEllipse(xyr, xy, aspectratio);
    setappdata(0,'EllipseFits',[xyr score]);
end

options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
    [x, score, exitflag, fitdetails] = fminsearch(@TryEllipse, xyr,options,  xy, aspectratio, trackfit);
[ score, rs] = TryEllipse(x, xy,aspectratio);
    x(5) = x(5);
ell = x;
    
function [cost, r] = TryEllipse(xyr, xy, aspectratio, trackfit)
%Tryellips uses radians/10 for the rotation, as this gives better
%convergence (!) 
%here idlist is true/false (1,0);
if nargin < 4
    trackfit = 0;
end
%aspectratio = 1;
xyr(4) = xyr(4)./aspectratio;
xy(:,2) = xy(:,2)./aspectratio;
[r, theta] = CalcClusterDistance(xyr,xy);
cost = sum((r-1).^2)./std(theta);
if trackfit
    X = getappdata(0,'EllipseFits');
    X(end+1,:) = [xyr cost];
     setappdata(0,'EllipseFits',X);
end



function [r, theta] = CalcClusterDistance(xyr, xy)
  xys = xyrotate(xy(:,1)-xyr(1),xy(:,2)-xyr(2),xyr(5));
  newxy = xys(:,1)./xyr(3) + i .* xys(:,2)./xyr(4);
  r = abs(newxy);
  theta = angle(newxy);
  
  
function ShowFits(xy, aspectratio)
    X = getappdata(0,'EllipseFits');
        
GetFigure('EllipseFits');
hold off;
plot(xy(:,1),xy(:,2),'.-');
hold on;

showfirst = 1;
if showfirst
    angles = 0:0.1:pi;
    for j = 1:length(angles)
        xyr([1 2]) = mean(xy);
        xyr(5) = angles(j); %angl
        xys = xyrotate(xy(:,1)-xyr(1),xy(:,2)-xyr(2),xyr(5));
        xyr([3 4]) = std(xys).*1.4;
        scores(j) = TryEllipse(xyr,  xy, aspectratio);
        xyrs(j,:) = xyr;
    end
    for j = 1:length(angles)
        DrawEllipse(xyrs(j,:));
    end
end
for j = 1:size(X,1)
    DrawEllipse(X(j,:));
end
