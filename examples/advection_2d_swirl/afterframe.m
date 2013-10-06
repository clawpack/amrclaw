axis([0 1 0 1]);
daspect([1 1 1]);

% yrbcolormap;

showpatchborders;

hidegridlines();
showpatchborders();
hidecontourlines;
setpatchborderprops(1:5,'linewidth',2);
showgridlines(5);

caxis([0 1]);

fprintf('%10s : %24.16f\n','qmin',qmin);
fprintf('%10s : %24.16f\n','qmax',qmax);

underover_colorbar;


clear afterframe;
