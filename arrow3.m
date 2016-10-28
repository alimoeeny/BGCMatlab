function arrow3(x,y,z,headangle,headlength,varargin)
% function arrow(x,y,headangle,headlength,varargin)
% Draws a single-headed skeleton arrow (not a filled-in triangle). 
% Input arguments are x,y (vectors as in plot command), headngle = angle between vector and the head of the arrow (degrees), headlength = length of arrow head
% May be followed by other paramter-value commands as in plot. 
% NB arrowhead looks good only if axes are equal.
%
% Example:
% arrow([3 4],[5 6],30,0.2,'col','r','linewidth',2)
%
%Modified from arrow.m by bgc 2016

hold_state = ishold;

plotcmds = varargin;
plot3(x,y,z,plotcmds{:});

hold on
arrowx = x(end)-x(1);
arrowy = y(end)-y(1);
arrowz = z(end)-x(1);
arrlength = sqrt(arrowx^2+arrowy^2+arrowz^2);
arrlength = sqrt(arrowx^2+arrowy^2);
cs = cos(headangle*pi/180);
sn = sin(headangle*pi/180);
% Far end
headx = (cs * arrowx + sn * arrowy) * headlength/arrlength;
heady = (-sn * arrowx + cs * arrowy) * headlength/arrlength;
plot3(x(end) - [0 headx], y(end) - [0 heady], z(end) + [0 0 ],plotcmds{:})
headx = (cs * arrowx - sn * arrowy) * headlength/arrlength;
heady = (sn * arrowx + cs * arrowy) * headlength/arrlength;
plot3(x(end) - [0 headx], y(end) - [0 heady], z(end) + [0 0 ],plotcmds{:})

if ~hold_state, hold off; end
