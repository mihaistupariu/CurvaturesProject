function [mat_NrN, mat_NV, mat_NT] = determine1RingNeigh(number_vertices, triang_mesh, vertices_ch)

%% Description
% Generates information on the 1-Ring neighborhood for each vertex on the
% basis of the triangle mesh topology
% Input.
%   number_vertices: the number of vertices
%   triang_mesh: the triangle meshes
%   vertices_ch: matrix indicating with an 1 the vertices on the convex hull
% Output.
%   mat_NrN: column matrix, indicating for each vertex the nr. of neighbors
%   mat_NV: row r of this matrix indicates the IDs of the vertices that are neighbors
%           to vertex r, by RESPECTING their order in the 1-R
%   mat_NT: row r of this matrix indicates the IDs of the triangles that are neighbors
%           to vertex r


%% Find the number of neighbours by using the 'nonordered' approach
mat_NrN=zeros(number_vertices,1);
for vv=1:number_vertices 
    % find the incident triangles
    neighbors=vertexAttachments (triang_mesh, vv);
    % transform the cell in array
    neighbors_mat=cell2mat (neighbors);
    % dim_neighbors is the number of adjacent vertices
    [~, dim_neighborsaux] = size(neighbors_mat); 
    if vertices_ch(vv,1)==0
        dim_neighbors=dim_neighborsaux; % if interior
    else 
        dim_neighbors=dim_neighborsaux+1; % if on the border
    end
    mat_NrN (vv,1)=dim_neighbors;
end

%% Size of matrices for neighbors and adjacent triangles
nr_max_vecini=max(mat_NrN);
mat_NV=zeros(number_vertices, nr_max_vecini);
mat_NT=zeros(number_vertices, nr_max_vecini);

%% Main loop
for vv=1:number_vertices
    [var1,var2]=determineOrderedNeighbors(triang_mesh, vertices_ch, vv);
    var_vf=mat_NrN(vv,1);
    if vertices_ch(vv,1)==0
        var_tr=var_vf;
    else 
        var_tr=var_vf-1;
    end
    mat_NV(vv,1:var_vf)=var1;
    mat_NT(vv,1:var_tr)=var2;
end



