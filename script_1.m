%% SCRIPT_1: 

%% Description 
% Handle triangles and provide triangle-relevant information for each 
% vertex for the computation of GC and MC (cf. Dyn; Meyer)
 
%% Initializations

% find the number of triangles: nr_triang
[nr_triang, b] = size (T); 

%% Triangle level
% Compute for the triangles the normal, the slope for each triangle in T 
T_normal=zeros(nr_triang,4);
% the area
T_area=zeros(nr_triang,1);
% the mixed area of triangle's vertices (Meyer et al., 2002)
T_mixedarea=zeros(nr_triang,3);
% the length of egdes 
T_edgeslength=zeros(nr_triang,3);
% the dot products
T_dots=zeros(nr_triang,3);
% the cosines 
T_cos=zeros(nr_triang,3);
% the angles
T_angles=zeros (nr_triang,3);


%% Vertex level
% Compute the vector of angular defects
% the number of adjacent triangles
V_adjacent_tri=zeros(nr_vf,1);
% the normal
V_normal_aux=zeros(nr_vf,3);
% the weighted normal (weights=areas)
V_weighted_normal_aux=zeros(nr_vf,3);
% the area
V_area=zeros(nr_vf,1);
% the mixed area contribution
V_mixedarea=zeros(nr_vf,1);
% the angle defect(first column in radians, second in deg)
V_AD_aux=zeros(nr_vf,1); 
% the vector sum cotan*vector (needed for mean curvature, MEYER)
V_mean_curvature_vector_aux=zeros(nr_vf,3);
% normal, weighted normal 
V_normal=zeros(nr_vf,3);
V_weighted_normal=zeros(nr_vf,3);

%% Main loop; all triangles are processed

  for vt = 1:nr_triang % take all triangles 
     % find the vertices 
     index1=T(vt,1); % first vertex
     index2=T(vt,2); %second vertex
     index3=T(vt,3); % third vertex
     M1=V(index1,1:3); % the first vertex, read from line index1 of matrix V
     M2=V(index2,1:3); % the second vertex, read from line index2 of matrix V
     M3=V(index3,1:3);   % the third point, read from line index3 of matrix V
     
     %% Information for triangles
     
     % the normal, the slope and the area of the triangle
      crossprod = cross (M2-M1, M3-M2); % compute the cross product
            if (crossprod(1,3)<0) % the last component: POSITIVE
                crossprod(1, 1:3)=-crossprod(1,1:3);
            end
        norma = norm (crossprod, 2);  % find the norm of the cross product
        crossprod(1,:)= crossprod (1,:)/norma; %normalize it
        cosslope=abs(crossprod(1,3));  % ABS of the third component is cos
    % the normal & the slope
    T_normal(vt, 1:3)=crossprod(1, 1:3);
    slope = acos (cosslope);
    T_normal(vt,4)= (slope*180)/pi;
    % the area of triangle #vt (half of norm of cross product)
    T_area(vt,1)= norma/2; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % length of the edges; by convention edge #i is opposite to vertex #i
    T_edgeslength (vt,1) = norm (M3-M2);
    T_edgeslength (vt,2) = norm (M1-M3);
    T_edgeslength (vt,3) = norm (M2-M1);
    % dot products; for vertex #i is on position #i
    T_dots (vt,1) = dot (M2-M1, M3-M1);
    T_dots (vt,2) = dot (M3-M2, M1-M2);
    T_dots (vt,3) = dot (M1-M3, M2-M3);
    % the cosines; for vertex #i is on position #i
    T_cos (vt,1)= T_dots(vt,1)/(T_edgeslength (vt,2)*T_edgeslength (vt,3));
    T_cos (vt,2)= T_dots(vt,2)/(T_edgeslength (vt,3)*T_edgeslength (vt,1));
    T_cos (vt,3)= T_dots(vt,3)/(T_edgeslength (vt,1)*T_edgeslength (vt,2));
    % the angles
    for i=1:3
        if abs(T_cos (vt,i))>= 1 
            disp ('ERROR1'); T_angles(vt,i)=0; 
        else        
        T_angles(vt, i)=acos (T_cos(vt,i));
        end
    end
    %  cotangents, needed later, too
    cotan1=T_cos(vt,1)/sqrt(1-T_cos(vt,1)^2);
    cotan2=T_cos(vt,2)/sqrt(1-T_cos(vt,2)^2);
    cotan3=T_cos(vt,3)/sqrt(1-T_cos(vt,3)^2);
    % the mixed area
    if ((T_cos(vt,1)>0&&T_cos(vt,2)>0)&& T_cos(vt,3)>0)
        T_mixedarea(vt,1)=(T_edgeslength(vt,2)^2*cotan2+T_edgeslength(vt,3)^2*cotan3)/8;
        T_mixedarea(vt,2)=(T_edgeslength(vt,3)^2*cotan3+T_edgeslength(vt,1)^2*cotan1)/8;
        T_mixedarea(vt,3)=(T_edgeslength(vt,1)^2*cotan1+T_edgeslength(vt,2)^2*cotan2)/8;
    end
    if T_cos (vt,1)<=0
        T_mixedarea(vt,1)=T_area(vt,1)/2;T_mixedarea(vt,2)=T_area(vt,1)/4;T_mixedarea(vt,3)=T_area(vt,1)/4;
    end
    if T_cos (vt,2)<=0
        T_mixedarea(vt,1)=T_area(vt,1)/4;T_mixedarea(vt,2)=T_area(vt,1)/2;T_mixedarea(vt,3)=T_area(vt,1)/4;
    end
    if T_cos (vt,3)<=0
        T_mixedarea(vt,1)=T_area(vt,1)/4;T_mixedarea(vt,2)=T_area(vt,1)/4;T_mixedarea(vt,3)=T_area(vt,1)/2;
    end
    
    %% Information transfer to vertices
    % number of adjacent triangles
    V_adjacent_tri(index1,1)=V_adjacent_tri(index1,1)+1; 
    V_adjacent_tri(index2,1)=V_adjacent_tri(index2,1)+1; 
    V_adjacent_tri(index3,1)=V_adjacent_tri(index3,1)+1;
    % the normal  
    V_normal_aux(index1,1:3)=V_normal_aux(index1,1:3)+T_normal(vt,1:3);
    V_normal_aux(index2,1:3)=V_normal_aux(index2,1:3)+T_normal(vt,1:3);
    V_normal_aux(index3,1:3)=V_normal_aux(index3,1:3)+T_normal(vt,1:3);
    % the weighted normal
    V_weighted_normal_aux(index1, 1:3)= V_weighted_normal_aux(index1, 1:3)+ T_area(vt,1)*T_normal(vt,1:3);
    V_weighted_normal_aux(index2, 1:3)= V_weighted_normal_aux(index2, 1:3)+ T_area(vt,1)*T_normal(vt,1:3);
    V_weighted_normal_aux(index3, 1:3)= V_weighted_normal_aux(index3, 1:3)+ T_area(vt,1)*T_normal(vt,1:3);        
    % the area
    V_area(index1,1)=V_area(index1,1)+T_area(vt,1); 
    V_area(index2,1)=V_area(index2,1)+T_area(vt,1); 
    V_area(index3,1)=V_area(index3,1)+T_area(vt,1); 
    % the mixed area
    V_mixedarea (index1,1)=V_mixedarea(index1,1)+T_mixedarea(vt,1);
    V_mixedarea (index2,1)=V_mixedarea(index2,1)+T_mixedarea(vt,2);
    V_mixedarea (index3,1)=V_mixedarea(index3,1)+T_mixedarea(vt,3);
    % the sum of the adjacent angles (in radians)
    V_AD_aux(index1, 1)= V_AD_aux(index1, 1)+T_angles (vt,1);
    V_AD_aux(index2, 1)= V_AD_aux(index2, 1)+T_angles (vt,2);
    V_AD_aux(index3, 1)= V_AD_aux(index3, 1)+T_angles (vt,3);
    % the product cotan*vector (needed for mean curvature /Meyer)
    V_mean_curvature_vector_aux(index1,:)= V_mean_curvature_vector_aux (index1,:)+cotan2*(M1-M3)+cotan3*(M1-M2);
    V_mean_curvature_vector_aux(index2,:)= V_mean_curvature_vector_aux (index2,:)+cotan3*(M2-M1)+cotan1*(M2-M3);
    V_mean_curvature_vector_aux(index3,:)= V_mean_curvature_vector_aux (index3,:)+cotan1*(M3-M2)+cotan2*(M3-M1);
  end  
      
%% Compute normal at each vertex
 
for ii=1:nr_vf
    % normal & weighted normal
    V_normal(ii,1:3)=V_normal_aux(ii,1:3)/(norm(V_normal_aux(ii,1:3)));
    V_weighted_normal(ii,1:3)= V_weighted_normal_aux(ii,1:3)/norm(V_weighted_normal_aux(ii,1:3));
end
 