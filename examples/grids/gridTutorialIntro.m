%% Overview of Grid Factory Routines
% In this tutorial, we will introduce you briefly to the grid structure and
% give an overview of some of the many grid factory routines that can be
% used to generate various types of grids.

%% Create a 2D grid
% We start by creating a simple 2D Cartesian grid. We then remove one cell
% and add geometry information to the resulting grid.
G = cartGrid([3,2]);
G = removeCells(G, 2);
G = computeGeometry(G);

% Plot the grid
newplot
plotGrid(G,'FaceColor',[0.95 0.95 0.95]); axis off;

%% Plot cell, face, and node numbers
% MRST uses a fully unstructured grid format in which cells, faces, and
% nodes, as well as their topological relationships, are represented
% explicitly. To illustrate, we extract the cell and face centroids as well as
% the coordinates of each node. These will be used for plotting the cells,
% faces and node indices, respectively.
c_cent = G.cells.centroids;
f_cent = G.faces.centroids;
coords = G.nodes.coords;

% Add circles around the centroids of each cell
hold on;
pltarg = {'MarkerSize',20,'LineWidth',2,'MarkerFaceColor',[.95 .95 .95]};

plot(c_cent(:,1), c_cent(:,2),'or',pltarg{:});

% Plot triangles around face centroids
plot(f_cent(:,1), f_cent(:,2),'sg',pltarg{:});

% Plot squares around nodes
plot(coords(:,1), coords(:,2),'db',pltarg{:});

legend({'Grid', 'Cell', 'Face', 'Node'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal')

% Plot cell/face centroids and nodes
txtargs = {'FontSize',12,'HorizontalAlignment','left'};
text(c_cent(:,1)-0.04, c_cent(:,2), num2str((1:G.cells.num)'),txtargs{:});
text(f_cent(:,1)-0.045, f_cent(:,2), num2str((1:G.faces.num)'),txtargs{:});
text(coords(:,1)-0.075, coords(:,2), num2str((1:G.nodes.num)'),txtargs{:});

title('Grid structure')
hold off;

%% Finding mappings between grid primitives
% The unstructured grid is built upon coordinates representing nodes, nodes
% representing faces and faces representing cells. To save memory, some
% attributes are run length encoded. For more information, see
% 'help grid_structure'.
% For instance, let us look up information for the faces.
faces =[ rldecode(1 : G.cells.num,diff(G.cells.facePos), 2).' G.cells.faces];
tag = {'West'; 'East'; 'South'; 'North'; 'Bottom'; 'Top'};
fprintf('Cell\tFace\tTag id\tTag\n');
for i=1:size(faces,1)
   fprintf(' %3d\t%3d\t%3d\t[%s]\n', faces(i,1:3), tag{faces(i,3)});
end

%% Neighborship is defined through faces
% G.faces.neighbors(i,:) contains the cells neighboring to face i. We
% define a new grid and plot the neighbors of face 10 in the new grid.
clf;
plotGrid(G,'FaceAlpha', 1, 'FaceColor', [0.95 0.95 0.95]); axis off;

f = 6;
plotGrid(G, G.faces.neighbors(f,:), 'FaceColor', 'Green')
text(f_cent(f,1)-0.1, f_cent(f,2), num2str(f),'FontSize',16,'Color','red');
% This also defines the boundary faces, since faces with only one neighbor
% is on the edge of the domain:
boundary = any(G.faces.neighbors==0,2);
facelist = 1:G.faces.num;
fprintf('Boundary faces: \n')
facelist( boundary) %#ok intentional display
fprintf('Internal faces: \n')
facelist(~boundary) %#ok intentional display

%% Generating synthethic grids using MRST
% There are many options for creating synthethic grids in MRST. For
% instance, the simple cartGrid already shown is a special case of
% tensorGrid:
G = tensorGrid((1:10).^3, 1:5);

clf;
plotGrid(G);

%% Triangulated grids
% We can generate a triangulated grid using existing triangulations, from
% MATLAB or elsewhere.

% Generate points
pts = rand(20,3).*repmat([10,10,1], 20, 1);
% Triangulate
T = delaunayn(pts);
G = tetrahedralGrid(pts, T);

clf, title('Tetrahedral grid')
plotGrid(G);
view(50,60), axis tight off

%% Triangular grids
% Generate a 2D triangle grid from the same data using the xy-coordinates
pts_2d = pts(:,1:2);
% Triangulate
T = delaunay(pts_2d(:,1), pts_2d(:,2));
G_tri = triangleGrid(pts_2d, T);

clf, title('Triangular grid')
plotGrid(G_tri); axis tight off

%% Extruded grid
% If we have an interesting 2D grid, we can easily extend it to 3D by using
% makeLayeredGrid. We will extend the previous triangle grid to 3 layers
G = makeLayeredGrid(G_tri, 3);
clf;
title('Extruded triangular grid')
plotGrid(G);
view(50,60), axis tight off

%% Explicit hexahedral grid
% Hexahedral grids can also be generated by lists of nodes and node
% indices. For valid node ordering, see help hexahedralGrid

H = [1  2  3  4  5  6  7  8; ... % Cell 1
    2  9 10  3  6 11 12  7]; ... % Cell 2


P = [  1  0  0.1860; ...
       1  1  0.1852; ...
       1  1  0.1926; ...
       1  0  0.1930; ...
       0  0  0.1854; ...
       0  1  0.1846; ...
       0  1  0.1923; ...
       0  0  0.1926; ...
       1  2  0.1844; ...
       1  2  0.1922; ...
       0  2  0.1837; ...
       0  2  0.1919]; ...

 G = hexahedralGrid(P, H);

 clf;
 plotGrid(G);
 axis tight
 view(40,40)

%% Grids can be manipulated after creation
% We can alter the attributes of a grid after creation. In this example we
% twist the grid slightly. One caveat: Any computed properties of the grid
% will not be altered. This means that computeGeometry must be called again
% after grid updates to get correct centroids, volumes, normals, etc.
G = cartGrid([10, 10]);
G_before = computeGeometry(G);
% Twist the coordinates to create a non-K-orthogonal grid.
G_after = twister(G);
G_after = computeGeometry(G_after);
clf;
plotGrid(G_after);
hold on
plot(G_before.cells.centroids(:,1), G_before.cells.centroids(:,2), 'Or')
plot(G_after.cells.centroids(:,1), G_after.cells.centroids(:,2), 'xb')
legend('Twisted grid', 'Original centroids', 'Modified centroids', ...
       'Location', 'NorthOutside', 'Orientation', 'horizontal')

%% Some grid routines produce corner-point input data
% Corner-point grids are usually represented on a special format given in
% ECLIPSE input files. MRST have a few routines that generate special grid
% models given on this input format, which in turn can be converted to MRST
% grids using processGRDECL, just as one would with a GRDECL file read
% using readGRDECL. For instance, here is a three layered structure which
% is easy to generate by creating pillars manually.

G = processGRDECL( threeLayers(10,10,5));
G = computeGeometry(G);
clf;
% Color the cells by the cell volume to show the layered structure.
plotCellData(G, G.cells.volumes,'EdgeColor','k');
view(120,10);
axis tight off

%% Grids generated by simpleGrdecl
% The function |simpleGrdecl| can generate different example grids based on
% optional parameters. The default parameter gives a simple, wavy grid.
grdecl = simpleGrdecl([20, 20, 5]);
G1 = processGRDECL(grdecl);
clf, plotGrid(G1), view(3), axis tight off;

%%
% By supplying a function handle, the grid becames faulted based on the
% values of the function at the fault line.
grdecl = simpleGrdecl([20, 20, 5], @(x) 0.05 * (sin(2*pi*x) - 1.5));
G2 = processGRDECL(grdecl);
clf, plotGrid(G2), view(3), axis tight off;

%%
% Create a flat grid with another fault function
grdecl = simpleGrdecl([20, 20, 5], @(x) 0.25*(x-0.5), 'flat', true);
G3 = processGRDECL(grdecl);
clf, plotGrid(G3), view(3), axis tight off;
%%
displayEndOfDemoMessage(mfilename)

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
