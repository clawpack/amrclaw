yrbcolormap;
daspect([1 1 1]);

showgridlines(1:2);

caxis([0.1 1]);

showpatchborders;
cv = linspace(0.2,0.9,8);
drawcontourlines(cv);

clear afterframe;
