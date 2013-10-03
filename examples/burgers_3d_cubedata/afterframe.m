yrbcolormap
caxis([0 1])
grid on
axis([0 1 0 1 0 1])
daspect([1 1 1])

showcubes;
setcubecolor('k',1);
setcubecolor('b',2);
setcubecolor('r',3);

colorbar;
camlight;
shg;

clear afterframe;
