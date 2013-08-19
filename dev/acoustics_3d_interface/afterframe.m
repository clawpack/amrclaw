yrbcolormap;
axis image;
grid on

showpatchborders;
showgridlines(1:2);

NoQuery = 0;
prt = false;
if (prt)
  fstr = framename(Frame,'frame0000','png','./');
  print('-dpng',fstr);
end;


clear afterframe;
