if PlotType ~= 4
   axis([0 1.2 0 0.5 0 0.5]);
   daspect([1 1 1]);
   end

if PlotType==1
  yrbcolormap
  end

if PlotType==5
  camlight left;
  grid off;
  end

clear afterframe
