daspect([1 1 1]);
colormap(jet);

showgridlines(1:2);

caxis([0.1 2]);

showpatchborders;
cv = linspace(0.2,1.6,16);
drawcontourlines(cv);

colorbar;

clear afterframe;
