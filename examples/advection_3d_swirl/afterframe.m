yrbcolormap
axis([0 1 0 1 0 1])
daspect([1 1 1]);

% hideslices;
% showslices('x',3);
% showslices('y',3);
% showslices('z',3);


showpatchborders;
setpatchbordercolor('k');
% setpatchborderprops(1,'linewidth',2,'color','k');  % new version only
% setpatchborderprops(2,'linewidth',2,'color','k');  % new version only
% setpatchborderprops(3,'linewidth',2,'color','k');  % new version only

showcubes;
setcubecolor('k',1);
setcubecolor('k',2);
setcubecolor('k',3);

showgridlines(1:2);

%cv = linspace(0,1,11);
%cv([1 end]) = [];
%drawcontourlines(cv);

h = surflight;

axis off;

shg;

clear afterframe;
