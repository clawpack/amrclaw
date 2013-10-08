axis([-1 1 -1 1]);
daspect([1 1 1]);

colormap(jet);

showpatchborders;
showgridlines(1:2);

cv = linspace(-1,1,21);
cv(cv == 0) = [];
drawcontourlines(cv);

caxis([-1 1]);

shg;

clear mapc2p
clear afterframe;
