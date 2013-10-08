axis([0 1 0 1]);
daspect([1 1 1]);

% Colorbar is set in 'underover';
% yrbcolormap;

showpatchborders;

hidegridlines();
showpatchborders();
hidecontourlines;
showgridlines(1:2);

caxis([0 1]);

fprintf('%10s : %24.16f\n','qmin',qmin);
fprintf('%10s : %24.16f\n','qmax',qmax);

qlo = 0;
qhi = 1;
under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
colorbar_underover(under_label,over_label);
 

shg;

clear underover;
clear afterframe;
