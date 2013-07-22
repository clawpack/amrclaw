axis([-1 1 -1 1 -1 1])
caxis([-.004 .004])
grid on

% plot dashed line showing where interface is:
hold on;
plot3([-1 1],[0 0],[0 0],'k--')
plot3([0 0],[-1 1],[0 0],'k--')
hold off

showpatchborders;
setpatchborderprops(1,'linewidth',2,'color','k');  % new version only
setpatchborderprops(2,'linewidth',2,'color','b');  % new version only
% setpatchbordercolor('r',1);
% setpatchbordercolor('b',2);

showcubes;
setcubecolor('k',1);
setcubecolor('b',2);


colorbar;

shg;

clear afterframe;
