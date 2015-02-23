if (PlotType == 1)
    % Assume presssure is being plotted.
    % rybcolormap;
    axis([0 2 0 2 0 2]);
    %caxis([0.2 0.6]);

    % for density:
    caxis([0.96 1.04]);

    daspect([1 1 1]);
    colorbar;
    
    view([136.5 28]);
    
    hideslices;
    showslices('x',1);
    showslices('y',1);
    showslices('z',1);
    
    showpatchborders;
    setpatchbordercolor('k');
    
    % showcubes;
    setcubecolor('k');
    
elseif (PlotType == 4)
   % Scatter plot   
   axis([0 2*sqrt(3) 0 qmax]);   
end


shg;
clear afterframe;
