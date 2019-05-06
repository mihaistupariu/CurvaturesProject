function [mat_2nrVf, mat_2vf, mat_2nrTr, mat_2tr, mat_2nrEdges, mat_2edges] =...
    determine2RingNeigh(vertices, triangulation, mat_nr, mat_vf, mat_tr, mat_edges)

%% Description
% Finds the 2-Ring Neighborhood starting from the 1R Neighborhood
% Input: 
%   vertices: matrix of vertices-->need their number;
%   mat_nr, mat_vf, mat_tr, mat_edges:  information about the 1-RN...
%   (nr of neighbours, adjacent vertices , adjacent triangles, adjacent
%   edges for each vertex)
% Output: 
%   information about the 2-RN (mat_2nrVf, mat_2vf, mat_2nrTr,
%   mat_2tr, mat_2nrEdges, mat_2edges); for each geometric element (vertices
%   / triangles edges) a matrix with their numbers followed by the elements -
%   see also below
 

%% Initializations
[nrVf,~]=size(vertices); % number of vertices
[nrTri,~]=size(triangulation); % number of triangles
[nrEdges,~]=size(mat_edges); %number of edges
mat_2nrVf=zeros(nrVf,1); % number of vertices in the 2-RN, per vertex
mat_2nrTr=zeros(nrVf,1); % number of triangles in the 2-RN, per vertex
mat_2nrEdges=zeros(nrVf,1); % number of edges in the 2-RN, per vertex
mat_2vf_aux=zeros(nrVf,nrVf); % vertices in the 2-RN, per vertex
mat_2tri_aux=zeros(nrVf,nrTri); % triangles in the 2-RN, per vertex
mat_2edges_aux=zeros(nrVf,nrEdges); % edges in the 2-RN, per vertex

%% Determine the vertices in the 2-Ring that are not in the 1-Ring
maxNrVec=0; % maximal number of neighbours
for vv=1:nrVf % take all vertices
    nrNeighv=mat_nr(vv,1); % the neighbours of v
    vect_aux=[mat_vf(vv,1:nrNeighv) vv]; % the triangles in the 1-RN
    neighs=[];
    % take all the neighbours
    for w=1:nrNeighv %take each neighbour
        ww=mat_vf(vv,w); % select the neighbour on position w!!!!!
        nrNeighw=mat_nr(ww,1); % the neighbours of v
        wect_aux=mat_vf(ww,1:nrNeighw); % the neighbours of w
        aux=setdiff(wect_aux, vect_aux); % the neighbours in the 2R induced by w
        neighs=union(neighs, aux); % add to the existing ones
    end
    % process the vertices found
           [~, nc]=size(neighs); % determine the size
    if nc==1 
        neighs=transpose(neighs); % transpose if column vector
    end
    [~, nc]=size(neighs);
    mat_2nrVf(vv,1)=nc;
    mat_2vf_aux(vv,1:nc)=neighs; 
    if nc>maxNrVec
        maxNrVec=nc; % if more than maximum remember this value
    end
end  
mat_2vf=mat_2vf_aux(1:nrVf,1:maxNrVec);

%% Determine the triangles in the 2-Ring that are not in the 1-Ring
maxNrVec=0; %maximal number of neighbours
for vv=1:nrVf % take all vertices
    nrNeighv=mat_nr(vv,1); % the neighbours of v
    vect_aux=mat_tr(vv,1:nrNeighv); % the triangles adjacent to v
    neighs=[];
    for w=1:nrNeighv % take all vertices in the 1-RN
        ww=mat_vf(vv,w); % select the neighbour on position w!!!!!
        wect_aux=cell2mat(vertexAttachments(triangulation, ww));
        aux=setdiff(wect_aux, vect_aux); % the triangles in the 2R induced by w
        neighs=union(neighs, aux); % add to the existing ones
    end
    % process the triangles found
    [~, nc]=size(neighs); % determine the size
    if nc==1 
        neighs=transpose(neighs); % transpose if column vector
    end
    [~, nc]=size(neighs);
    mat_2nrTr(vv,1)=nc;
    mat_2tri_aux(vv,1:nc)=neighs; 
    if nc>maxNrVec
        maxNrVec=nc; % if more than maximum remember this value
    end
end
mat_2tr=mat_2tri_aux(1:nrVf,1:maxNrVec);

%% Determine the edges in the 2_ring that are not in the 1-Ring
maxNrVec=0; %maximal number of neighbours
for vv=1:nrVf
    nrNeighv=mat_nr(vv,1); % the neighbours of v
    vect_aux=mat_edges(vv,1:nrNeighv); % the edges adjacent to v
    neighs=[];
    for w=1:nrNeighv % take all vertices in the 1-RN
        ww=mat_vf(vv,w); % select the neighbour on position w!!!!!
        nrNeighw=mat_nr(ww,1); % the neighbours of v
        wect_aux=mat_edges(ww,1:nrNeighw);
        aux=setdiff(wect_aux, vect_aux); % the triangles in the 2R induced by w
        neighs=union(neighs, aux); % add to the existing ones
    end
    % process the edges found
           [~, nc]=size(neighs); % determine the size
    if nc==1 
        neighs=transpose(neighs); % transpose if column vector
    end
    [~, nc]=size(neighs);
    mat_2nrEdges(vv,1)=nc;
    mat_2edges_aux(vv,1:nc)=neighs; 
    if nc>maxNrVec
        maxNrVec=nc; % if more than maximum remember this value
    end
end
mat_2edges=mat_2edges_aux(1:nrVf,1:maxNrVec);
end

