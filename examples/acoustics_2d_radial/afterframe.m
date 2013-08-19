if (PlotType == 1)
  daspect([1 1 1]);
  rybcolormap;
  showgridlines(1:2);
  showpatchborders;

  caxis([-2 2]);
  colorbar;
elseif (PlotType == 4)
  hold on;
  dir = './1drad/';
  dim = 1;
  [amrdata1d,t1d] = readamrdata(dim,Frame,dir);
  % [q1d,x1d] = plotframe1ez(amrdata1d,mq,'b-');

  ph = getlegendinfo;
  lh = legend(ph,{'level 1','level 2','level 3'});
  set(lh,'fontsize',16);

  hold off;
end


prt = false;
if (prt)
  fname = framename(Frame,'radial0000','png');
  print('-dpng',fname);
end



clear afterframe;
