yrbcolormap
grid on
axis([0 1 0 1 0 1])
daspect([1 1 1]);



showpatchborders;
setpatchborderprops(1,'linewidth',2,'color','k');  % new version only
setpatchborderprops(2,'linewidth',2,'color','b');  % new version only
setpatchborderprops(3,'linewidth',2,'color','r');  % new version only
%setpatchbordercolor('k',1);
setpatchbordercolor('b',2);
setpatchbordercolor('r',3);


%showcubes;
setcubecolor('k',1);
setcubecolor('b',2);
setcubecolor('r',3);


shg;
clear afterframe;
