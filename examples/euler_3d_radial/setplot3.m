
%  setplot3.m
%  called in plotclaw3 before plotting to set various parameters

OutputDir = '_output';

setopengl;
setviews  % set viewpoints so that view(xSlice), for example, can be used.

% PlotType = 4;                % type of plot to produce:
PlotType = 1;                % type of plot to produce:
			     % 1 = pcolor on slices (with optional contours)
			     % 2 = contour lines in 3d on transparent slices
			     % 3 = Schlieren plot on slices
			     % 4 = scatter plot of q vs. r
			     % 5 = isosurface plot (at given levels)

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = 'pressure';      % name of m-file mapping data to q
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 6;               % max number of AMR levels

PlotData =  [1 1 1 0 0 0];       % Data on refinement level k is plotted only
			         % if k'th component is nonzero
PlotGrid =  [1 1 0 0 0 0];       % Plot grid lines on each level?
PlotGridEdges =  [1 0 0 0 0 0];  % Plot edges of patches of each grid at
                                 % this level on slices?
PlotCubeEdges = [0 0 0 0 0 0];   % Plot edges of cube of refinement patch at
                                 % this level?


% ContourValues is a vector of contour lines that can be used with
% PlotType = 1,2,3.   Empty ==> no contour lines drawn
  ContourValues = [];

% The next three parameters are vectors of x,y,z coordinates of 2d slices
% to be displayed for PlotType = 1,2,3.
s = linspace(0,2,6);
  xSliceCoords = s;
  ySliceCoords = s;
  zSliceCoords = s;

% For PlotType = 4 (Scatter plot)
% Use default x0, y0, z0;
x0 = 0;
y0 = 0;
z0 = 0;
ScatterStyle = setplotstyle('b.','g+','ko');

% For Isosurfaces (which will plotted automatically for PlotType <= 3)
% Note:  Lengths of SurfTransparency and SurfColors must greater than or
%        equal to the length of SurfValues..

IsosurfValues = [];
