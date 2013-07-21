if PlotType == 1
  % Assume presssure is being plotted.
  % rybcolormap;
  axis([0 2 0 2 0 2]);
  caxis([0.2 0.6]);
  daspect([1 1 1]);
  colorbar;
end

if PlotType==4
   % scatter plot.  Also plot rad1d solution for comparison
   [rad1d,tref] = readamrdata(1,Frame,'../rad1d/');
   if (abs(tref  - t) > 1e-8)
     error('Times for reference solution and 3d solution not compatible');
   end
   hold on;
   if (UserVariable == 1)
     userfile = 'pressure';
   else
     userfile = '';
   end;
   set(gca,'Box','on');
   set(gca,'YLimMode','auto');
   [qref,xref,p] = plotframe1ez(rad1d,mq,'r-',userfile);
   set(p,'LineWidth',2.5);
   [h_amr,str_amr] = getlegendinfo;
   legend([h_amr,p],{'Level 1 (mx = 20)','Level 2 (mx = 40)',...
       'Level 3 (mx = 80)','1d solution (mx = 500)'});
 end

 clear afterframe;
