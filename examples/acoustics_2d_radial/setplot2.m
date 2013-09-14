OutputDir = './_output';

PlotType = 4;                % type of plot to produce:
% 			     - 1 = pcolor on slices (with optional contours)
% 			     - 2 = contour lines in 2d on white slices
% 			     - 3 = Schlieren plot on slices
% 			     - 4 = scatter plot of q vs. r

mq = 1;                      % which component of q to plot
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
Manifold = 0;
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 6;               % increase if using amrclaw with more levels
PlotData =  [1 1 1 1 1 1];   % Data on refinement level k is plotted only if
			     % k'th component is nonzero
PlotGrid =  [0 0 0 0 0 0];   % Plot grid lines on each level?

PlotGridEdges =  [1 1 1 1 1 1];  % Plot edges of patches of each grid at
                                 % this level?


ContourValues = [];  % Set in afterframe using drawcontourlines

%---------------------------------

UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = ' ';      % name of m-file mapping data to q

% Symbols for scatter plot.  Set a different symbol for each level
ScatterStyle = setplotstyle('bo','r*','gx');
x0 = 0;
y0 = 0;
