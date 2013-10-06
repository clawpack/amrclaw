function uo = underover()

% UNDEROVER returns structure for assigning colors to under/overshoots
%
% UO = UNDEROVER().  The user needs to define the fields as follows
%
%          color_under     the rgb triple used to color undershoots
%          color_over      the rgb triple used to color overshoots
%          value_lower     lower value of exact solution
%          value_upper     upper value of exact solution
%          tol             tolerance for under/over shoots
%          colormap        colormap for values in [qlow-tol,qhi+tol];
%
% Values are considered in the overshoot/undershoot region if they satisfy
%
%                 q < value_lower - tol
% or
%                 q > value_upper + tol
%
% Values not in the under/overshoot region are mapped into the colormap
% according to
%
%                    slope = (q-qlow)/(qhi-qlow)
%
%                    idx = 1 + slope*(nmax-1)
%
% An extended colormap is created as :
%
%                extended_colormap = [color_under; colormap; color_over];
%
% Using the extended colormap, the 'cdata' property of the finite volume
% patch is set to an rgb triple.  Values in the overshoot region are set to
% color_over, values in the undershoot region are set to color_under, values
% in [value_lower-tol, value_lower] are set to the first color in
% 'colormap',  values in [value_upper, value_upper+tol] are set to the last
% values in 'colormap', and values in [value_lower, value_upper] are
% linearly scaled into the colormap.
%
% This function is called from SETCOLORS
%
% See also SETCOLORS
%

yrbcolormap;
cm = colormap(jet(64));

uo = struct('color_under',[0 1 1],...
            'color_over',[1 0 1],...
            'value_lower',0, ...
	    'value_upper',1,...
	    'tol',1e-4,...
	    'colormap',cm);
