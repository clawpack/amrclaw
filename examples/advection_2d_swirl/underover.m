function uo = underover()
%UNDEROVER returns structure for assigning colors to under/overshoots
% 
%   For complete help, type 
% 
%   >> help visclaw/src/matlab/underover
%

yrbcolormap;
cm = colormap;

uo = struct('color_under',[0 1 1],...
    'color_over',[1 0 1],...
    'value_lower',0, ...
    'value_upper',1,...
    'tol',1e-6,...
    'colormap',cm);
