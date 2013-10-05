yrbcolormap
axis([0 1 0 1 0 1])
daspect([1 1 1]);

hideslices;
showslices('x',6);
showslices('y',6);
showslices('z',6);


showpatchborders;
setpatchborderprops(1,'linewidth',2,'color','k');  % new version only
setpatchborderprops(2,'linewidth',2,'color','b');  % new version only
setpatchborderprops(3,'linewidth',2,'color','r');  % new version only
%setpatchbordercolor('k',1);
setpatchbordercolor('k',2);
setpatchbordercolor('k',3);


showcubes;
setcubecolor('r',1);
setcubecolor('b',2);
setcubecolor('k',3);
hidecubes(1:2);

showgridlines(2:3);

cv = linspace(0,1,11);
cv([1 end]) = [];
drawcontourlines(cv);

axis off;

shg;

clear afterframe;
