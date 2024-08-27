close all
clear all
clc

% Mesh Class
mesh = MeshClass();
% number of elements
mesh.Nel = uint16([5 5 2]);
% domain size
geometry.Dimension = [1 1 0.25];
mesh.Mesher(geometry);

%xmin xmax; yminy max; zmin zmax
region = [0.8 1.2; 0.39 0.7;0 1.2];
%refine
mesh.RefMesh(region);

%plot
mesh.PlotMesh(figure)
axis square
% view(-27,17)
view(2)
%export to gmsh
mesh.MeshWrite('export.msh')
% it should have prod mesh.Nel
numel(unique(mesh.Father))
% all elements
numel(mesh.Father)
%number of created elements in the refinemnt
numel(mesh.Father)-numel(unique(mesh.Father))

