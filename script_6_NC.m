%% SCRIPT_6_NC: 

%% Description
% Implementation of the Normal Cycle approach(Cohen-Steiner & Morvan, 2003; Alliez et al, 2003) 
% Computes the 3x3 tensor from Cohen-Steiner & Morvan; Alliez et al
 
%% Initialization
tensorNC_edges=zeros(3,3, nr_edges); 
tensorNC_vertices_1R=zeros(3,3, nr_vf); 
tensorNC_vertices_2R_partial=zeros(3,3, nr_vf); 
[~,V_1R_Neigh_EdgesNr]=size(V_1R_Neigh_Edges);
[~,V_1R_Neigh_TriNr]=size(V_1R_Neigh_Tri);
% Curvatures
V_GC_NC_1R=zeros(nr_vf, 1);  
V_MC_NC_1R=zeros(nr_vf, 1);  
V_GC_NC_2R=zeros(nr_vf, 1);  
V_MC_NC_2R=zeros(nr_vf, 1);  

%% First loop: the edge loop; for each edge construct the 3x3 matrix 
for ee=1:nr_edges
        % find the vertices adjacent to the edge
        i1=E(ee,1);
        i2=E(ee,2);
        
        % %%%%% Edge measures (direction, length, versor, tensor)
        % compute characteristic measures of the current edge
        edge_dir=V(i1,1:3)-V(i2,1:3); % the direction of the edge
        edge_length=norm(edge_dir); % its length
        edge_versor=edge_dir/edge_length;
        edge_tensor=transpose(edge_versor)*edge_versor;
        
        % %%%%% Adjacent triangles (signed angle) 
        %find the adjacent triangles
        adjacent_tri=edgeAttachments (T, i1,i2);
        % transform the cell in array; this is the matrix of IDs of adj. tr
        adjacent_mat=cell2mat (adjacent_tri);
        % dim_neighbors is the number of adjacent triangles
        [~, dim_neighbors] = size(adjacent_mat);
        % if the edge has two adjacent triangles, compute the normal deviation and the contribution; else it is zero
        if (dim_neighbors==2)
            %find the adjacent triangles as j1, j2 and compute the angle of normals
            j1=adjacent_mat(1,1); vect_aux1=T_normal(j1, 1:3); % the first triangle and its normal
            j2=adjacent_mat(1,2); vect_aux2=T_normal(j2, 1:3); % the second triangle and its normal
            % indices of the vertices adjacent to the edge E
            i3=setdiff(T(j1,:), [i1 i2]); 
            i4=setdiff(T(j2,:), [i1 i2]);
            % get coordinates of the corresponding vertices
            varf1=V(i1,:); varf2=V(i2,:); varf3=V(i3,:); varf4=V(i4,:);
            % sign of the edge
            signofEdge=determineEdgeSign(varf1, varf2, varf3, varf4);
            dv1dv2=0;
            dv1v2=dot(vect_aux1, vect_aux2);
            % compute the normal deviation 
            if abs(dv1v2)<=1 % this should be a cosine, since both normals are true normals
                edge_norm_dev=acos (dv1v2); 
            else
                    if dv1dv2>1
                        edge_norm_dev=0;
                    else
                        edge_norm_dev=acos(-1);
                    end 
            end
            % the signed angle
            signed_angle=signofEdge*edge_norm_dev;
        else  
            % for other edges signed angle is 0
            signed_angle=0;
        end
        % %%%%%% MAIN COMPUTATION: Formula Cohen-Steiner&Morvan; Alliez et al
        tens_aux=signed_angle*edge_length*edge_tensor;
        tensorNC_edges(:,:,ee)=tens_aux;
end


%% Second loop: determine tensors for all vertices

% %%%%%% For the 1-Ring Neighborhood
for vv=1:nr_vf % take all vertices
     for ee_crt=1:V_1R_Neigh_EdgesNr % take all edge indices
         if V_1R_Neigh_Edges(vv,ee_crt)~=0 % not zero--meaningful
            edge_crt=V_1R_Neigh_Edges(vv,ee_crt); % find the edge corresponding to the position ee_crt
            tensorNC_vertices_1R(:,:,vv)=tensorNC_vertices_1R(:,:,vv)+tensorNC_edges(:,:,edge_crt);  % add the contribution of edge_crt to the tensor in vv
         end
     end
     % compute surrounding area
     area_ngh=0; % initialize neighbourhing area
     for tt_crt=1:V_1R_Neigh_TriNr % take all triangle indices
         if V_1R_Neigh_Tri(vv,tt_crt)~=0 % not zero--meaningful
             tri_crt=V_1R_Neigh_Tri(vv,tt_crt); % find the corresponding triangle
             ind1=T(tri_crt,1); ind2=T(tri_crt,2); ind3=T(tri_crt,3);  % indices of the vertices for the triangle tri_crt
             P1=V(ind1, 1:3); P2=V(ind2, 1:3); P3=V(ind3, 1:3); % coordinates of the vertices
             [area_tri_crt,~]=determineMeasuresTriangle(P1,P2,P3); %area 
             area_ngh=area_ngh+area_tri_crt;
         end
     end
     tensorNC_vertices_1R(:,:,vv)=tensorNC_vertices_1R(:,:,vv)/area_ngh; 
end

% %%%%%% For the 2-Ring Neighborhood
for vv=1:nr_vf % take all vertices
     for ee_crt=1:V_2R_Neigh_EdgesNr % take all edge indices for the 2-RN
         if V_2R_Neigh_Edges(vv,ee_crt)~=0 % not zero--meaningful
            edge_crt=V_2R_Neigh_Edges(vv,ee_crt); % find the edge corresponding to the position ee_crt
            tensorNC_vertices_2R_partial(:,:,vv)=tensorNC_vertices_2R_partial(:,:,vv)+tensorNC_edges(:,:,edge_crt);  % add the contribution of edge_crt to the tensor in vv
         end
     end
     % compute surrounding area
     area_ngh=0; % initialize neighbourhing area
     for tt_crt=1:V_2R_Neigh_TriNr % take all triangle indices
         if V_2R_Neigh_Tri(vv,tt_crt)~=0 % not zero--meaningful
             tri_crt=V_2R_Neigh_Tri(vv,tt_crt); % find the corresponding triangle
             ind1=T(tri_crt,1); ind2=T(tri_crt,2); ind3=T(tri_crt,3);  % indices of the vertices for the triangle tri_crt
             P1=V(ind1, 1:3); P2=V(ind2, 1:3); P3=V(ind3, 1:3); % coordinates of the vertices
             [area_tri_crt,~]=determineMeasuresTriangle(P1,P2,P3); %area 
             area_ngh=area_ngh+area_tri_crt;
         end
     end
     tensorNC_vertices_2R_partial(:,:,vv)=tensorNC_vertices_2R_partial(:,:,vv)/area_ngh; 
end
tensorNC_vertices_2R=tensorNC_vertices_1R+tensorNC_vertices_2R_partial;


%% Third loop: extract eigenvalues and compute 
 
% %%%% For the 1-Ring Neighborhood
for vv=1:nr_vf
    eigenv=zeros(2,1);
    tensor_aux=tensorNC_vertices_1R(:,:,vv);
    eig_aux=eig(tensor_aux);
    [minAbsEig,~]=find (abs(eig_aux)==min(abs(eig_aux)));
    indiceMinAbsEig=minAbsEig(1,1);
      indCrtRun=0;
      for indCrt=1:3
          if (indCrt~=indiceMinAbsEig)
              indCrtRun=indCrtRun+1;
              eigenv(indCrtRun,1)=eig_aux(indCrt,1);
          end
      end 
      V_GC_NC_1R(vv,1)=eigenv(1,1)*eigenv(2,1);
      V_MC_NC_1R(vv,1)=-(eigenv(1,1)+eigenv(2,1))/2;
end

% %%%% For the 2-Ring Neighborhood
for vv=1:nr_vf
    eigenv=zeros(2,1);
    tensor_aux=tensorNC_vertices_2R(:,:,vv);
    eig_aux=eig(tensor_aux);
    [minAbsEig,~]=find (abs(eig_aux)==min(abs(eig_aux)));
    indiceMinAbsEig=minAbsEig(1,1);
      indCrtRun=0;
      for indCrt=1:3
          if (indCrt~=indiceMinAbsEig)
              indCrtRun=indCrtRun+1;
              eigenv(indCrtRun,1)=eig_aux(indCrt,1);
          end
      end 
      V_GC_NC_2R(vv,1)=eigenv(1,1)*eigenv(2,1);
      V_MC_NC_2R(vv,1)=-(eigenv(1,1)+eigenv(2,1))/2;
end

