axis([-1 1 -1 1 -1 1])
caxis([-.004 .004])
grid on

% plot dashed line showing where interface is:
hold on;
plot3([-1 1],[0 0],[0 0],'k--')
plot3([0 0],[-1 1],[0 0],'k--')
hold off
