%% SCRIPT_2:

%% Description
% Handle edges and provide edge-relevant information for
% each vertex for the computation of GC and MC (Dyn; Meyer)

 

 
%% Initializations
% the 'edges contribution'
V_edge_contrib=zeros(nr_vf,1);
% process all edges
 

%% Main loop; all edges are processed 

for current_edge=1:nr_edges 
    %% Variation along the edge
    % find the vertices adjacent to the edge
        i1=E(current_edge,1);
        i2=E(current_edge,2);
    % compute characteristic measures of the current edge
        edge_dir=V(i1,1:3)-V(i2,1:3); % the direction of the edge
        edge_length=norm(edge_dir); % its length
    %find the adjacent triangles
        adjacent_tri=edgeAttachments (T, i1,i2);
        % transform the cell in array; this is the matrix of IDs of adj. tr
        adjacent_mat=cell2mat (adjacent_tri);
        % dim_neighbors is the number of adjacent triangles
        [dim_neighbors1, dim_neighbors] = size(adjacent_mat);
 
    %  lung(edge)*variation of normal*sign of edge
    % if the edge has two adjacent triangles, compute the normal deviation
    % and the contribution; else it is zero
        if (dim_neighbors==2)
            %find the adjacent triangles as j1, j2 and compute the angle of
            %normals
            j1=adjacent_mat(1,1); vect_aux1=T_normal(j1, 1:3); % the first triangle and its normal
            j2=adjacent_mat(1,2); vect_aux2=T_normal(j2, 1:3); % the second triangle and its normal
            % Introduce the sign and convex/concave edges 
            % indices of the vertices adjacent to the edge E
            i3=setdiff(T(j1,:), [i1 i2]); 
            i4=setdiff(T(j2,:), [i1 i2]);
            % get coordinates of the corresponding vertices
            varf1=V(i1,:); varf2=V(i2,:); varf3=V(i3,:); varf4=V(i4,:);
            % sign of the edge
            signofEdge=determineEdgeSign(varf1, varf2, varf3, varf4);
            dv1v2=dot(vect_aux1, vect_aux2);
            % compute the normal deviation 
            if abs(dv1v2)>1 % this should be a cosine, since both normals are true normals
                edge_norm_dev=0;
            else
            edge_norm_dev=acos (dv1v2);
            end
            % add the contribution of the edge to eaxh vertex; the sign of
            % the edge is also included
            V_edge_contrib(i1,1)=V_edge_contrib(i1,1)+(edge_length*edge_norm_dev*signofEdge)/4;
            V_edge_contrib(i2,1)=V_edge_contrib(i2,1)+(edge_length*edge_norm_dev*signofEdge)/4;   
        end
 
   
end

 
