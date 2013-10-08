colormap(jet);
axis([-1 1 -1 1 -1 1])
caxis([-0.2 .5])
grid on

% plot dashed line showing where interface is:
hold on;
plot3([-1 1],[0 0],[0 0],'y','linewidth',3)
plot3([0 0],[-1 1],[0 0],'y','linewidth',3)
hold off

showgridlines(1:2);

showpatchborders;
setpatchbordercolor('k');

cv = linspace(-0.005,0.005,20);
drawcontourlines(cv);
setcontourlineprops('color','w');

colorbar;

shg;

clear afterframe;
