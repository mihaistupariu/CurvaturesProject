%% SCRIPT_0: 

%% Description 
% Preprocess the input point cloud: generate the triangle mesh and get
% topological data
% The point records can be found in the matrix V. 
% Columns of V: 1-3 - coordinates
 
 
%% Initialization
[nr_vf,~]=size(V);
V_convexhull=zeros(nr_vf,1);

%% Generate the triangle mesh / convex hull
% Generate the Delaunay triangulation
T = delaunayTriangulation (V(:,1), V(:,2));  

% Determine the convex hull of T and associated values
convex_hull = convexHull(T);
[dch1, dch2]=size (convex_hull);

% Determine the size of the border and number of "interior" vertices (relevant for statistics)
size_border=dch1-1;
nr_vf_valid=nr_vf-size_border; 
% Determine the vertices arising on boundary of the convex hull (0=NO; 1=YES)
ii =1 ; 
while ii < dch1 
    xx = convex_hull(ii,1);
    V_convexhull(xx,1)=1;
    ii = ii+1;
end; 

%% Determine the 1-Ring Neighborhood for each vertex
%  Transfer topological information related to vertices / triangles to matrices
[V_1R_NeighNr, V_1R_Neigh, V_1R_Neigh_Tri]=determine1RingNeigh (nr_vf, T, V_convexhull);

% Information about edges
E=edges(T); %generate the edges: E is the matrix of edges
[nr_edges,~]=size(E); %find the size of E
[~,max_nr_neigh]=size(V_1R_Neigh); % maximum number of neighbours
V_1R_Neigh_Edges=zeros(nr_vf,max_nr_neigh); % the matrix of edge neighbours
for vv=1:nr_vf
    [nrNeighEdge, neighEdge]=determineAdjacentEdges(E,vv);
    V_1R_Neigh_Edges(vv,1:nrNeighEdge)=neighEdge;
end
 
%% Determine the 2-Ring Neighborhood for each vertex
[V_2R_NeighNr, V_2R_Neigh,V_2R_Neigh_TriNr, V_2R_Neigh_Tri,  V_2R_Neigh_EdgesNr, V_2R_Neigh_Edges]=...
    determine2RingNeigh(V,T,V_1R_NeighNr, V_1R_Neigh,V_1R_Neigh_Tri,V_1R_Neigh_Edges);