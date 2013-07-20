yrbcolormap
grid on
axis([0 1 0 1 0 1])

if (Frame == 0)
  camlight;
end;

% Reduce size of surface mesh;
reducesurf(0.6);
showsurfmesh;

clear afterframe;
