% add directory where Hybmesh.m is located
addpath('../../../../build/bindings/m');

% set hybmesh executable
Hybmesh.hybmesh_exec_path('../../../../src/py');

% set directory where shared library is located
Hybmesh.hybmesh_lib_path('../../../../build/bin/');

% initialize hybmesh object
hm = Hybmesh();

% set verbose stdout output
hm.stdout_verbosity(3);

% define set of points bounding triangulated area
% first equals last for closed contours.
points = [0.0, 0.0;
          0.5, 0.1;
          1.0, 0.2;
          1.3, 0.7;
          1.1, 1.3;
          0.3, 1.4;
          0.0, 1.0;
          0.0, 0.0];

% create a bounding contour
cont = hm.create_contour(points);

% perform contour parition using reference points:
% at point closest to [0, 0] recommended size is 0.05,
% at point closest to [1, 1] recommended size is 0.2.
domain = hm.partition_contour_ref_points(cont, [0.05, 0.2], [0, 0; 1, 1]);

% add two inner points which should present in result
inner_points = [0.5, 0.5; 0.7, 0.7];

% triangulate domain with two additional inner points
% with recommended step sizes 0.01, 0.03 respectively.
grid = hm.triangulate_domain(domain, [], [0.01, 0.03], inner_points);

% get grid vertices as plain [x0, y0, x1, y1, ....] array
pts = grid.raw_vertices();

% get cell-vertices table as plain array of vertex indicies where
% each three represent a triangle. Indexation starts with zero.
tris = grid.raw_tab('cell_vert');

% modify grid data to be able to plot it using triplot function:
% 1. slice pts array into separate x and y arrays
x=pts(1:2:end);
y=pts(2:2:end);
% 2. reshape tris into 2D matrix and add unity to fit native indexation.
tris = reshape(tris, 3, [])' .+ 1;

% save grid to vtk
hm.export_grid_vtk(grid, "out2.vtk");

% explicitly free hybmesh handle
hm.delete();

% plot data
triplot(tris, x, y);

% Enter keyboard interactive mode.
% This stops script from quiting.
keyboard;
