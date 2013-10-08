
OutputDir = '_output';

PlotType = 1;                % type of plot to produce:
			     % 1 = pseudo-color (pcolor)
			     % 2 = contour
			     % 3 = Schlieren
			     % 4 = scatter plot of q vs. r

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = ' ';      % name of m-file mapping data to q
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
Manifold = 0;
MaxFrames = 1000;            % max number of frames to loop over
ReadBlockNumber = 0;
MaxLevels = 6;
PlotData =  [1 1 1 1 1 1];   % Data on refinement level k is plotted only if
			     % k'th component is nonzero
PlotGrid =  [1 1 1 1 1 0];   % Plot grid lines on each level?

PlotGridEdges =  [1 1 1 1 1 1];  % Plot edges of patches of each grid at
                                 % this level?

% Set to either a scalar, for automatic contours or a vector of contour levels.
ContourValues = [];

% Show undershoots and over shoots.
ShowUnderOverShoots = 1;


%---------------------------------

% for scatter plot (PlotType==4):
% plot q(i,j) vs. r(i,j) = (x(i,j)-x0)^2 + (y(i,j)-y0)^2
  x0 = 0;
  y0 = 0;
  ScatterStyle = setplotstyle('bx','rx','gx');
