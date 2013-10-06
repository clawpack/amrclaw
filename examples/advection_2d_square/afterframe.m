daspect([1 1 1]);
yrbcolormap;

cv = linspace(0,1,6);
cv([1 end]) = [];
drawcontourlines(cv);

showgridlines(1:2);
showpatchborders;

shg

clear afterframe;
