function [mat_IDs, neighbors_mat] = determineOrderedNeighbors(triangulation, vertices_convhull, index_vf)
%% Description
% The function determines the 1-ring neighbours of a FIXED vertex in a (Delaunay)
% triangulation. The list of the vertices reflects their adjacencies in the Delaunay triangulation. 
% Input.
%   triangulation: the triangulation considered
%   vertices_convhull: matrix indicating with an 1 the vertices on the convex hull
%   index_vf: index of the considered vertex for which the neighbors are determined 
% Output.
%   mat_IDs: matrix of IDs of the neighbours in the 1-ring ngbh 
%	neighbors_mat: matrix of the neighboring triangles
 

%% Function that correctly 'links' the vertices of a new triangle (keeps the
% order correct)
function [vect_final]=unionLink (vect_init, vect_new, indice)
    [naj, ~]=size(vect_init);
    % remove from vect_new indice, which is a common vertex
    vect_aux2=setdiff(vect_new,indice);
    if (naj > 1)
        % the last two elements of vect_init
        vect_aux1(1,1)= vect_init(naj-1,1);
        vect_aux1(2,1)= vect_init(naj,1);
        % the first (naj-2) positions remain unchanged
        vect_final(1:naj-2,1)=vect_init(1:naj-2,1);
        % position (naj-1) is the difference between init and new            % new 
        vect_final(naj-1,1)= setdiff (vect_aux1, vect_aux2);
        % position naj is the interesection (common vertex); it links
        vect_final(naj,1)= intersect(vect_aux1, vect_aux2); 
        % position (naj+1) is the difference between new and init
        vect_final(naj+1,1)= setdiff (vect_aux2, vect_aux1);
    else 
        vect_final=vect_aux2;
    end
end

%% Find the incident triangles
neighbors=vertexAttachments (triangulation, index_vf);
% transform the cell in array; this is the matrix of IDs of adj. tr
neighbors_mat=cell2mat (neighbors);

%% dim_neighbors is the number of adjacent triangles
  [~, dim_neighborsaux] = size(neighbors_mat); 
    if vertices_convhull(index_vf)==0
        dim_neighbors=dim_neighborsaux; % if interior
    else 
        dim_neighbors=dim_neighborsaux+1; % if on the border
    end
    
%% Determine the matrix of neighbors
    mathelp =[];
    for vartr=1:dim_neighborsaux %take all adjacent triangles
        nr_tr=neighbors_mat(1,vartr); % find the ID of the triangle
        mathelp=unionLink (mathelp, transpose(triangulation(nr_tr,:)), index_vf); % attach the IDs of neighbors to the list
    end
   [dimhelp,~]=size(mathelp);
   if vertices_convhull (index_vf,1)==0
        mat_IDs=mathelp(1:dimhelp-1,:); % remove the last vertex for interior vertices
   else 
        mat_IDs=mathelp(1:dimhelp,:); %keep all vertices for border
   end
   [dim_neighbors_true,~]=size (mat_IDs);
    if (dim_neighbors~=dim_neighbors_true)
        disp ('ERROR NR NEIGHBOURS');
        disp (index_vf);
    end
 
    
end

