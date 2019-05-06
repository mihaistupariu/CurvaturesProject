%% SCRIPT_4_ET_TA: 
% Handle vertices and compute
% GC and MC for Euler's Theorem ET (Watanabe & Belyaev, 2001)
% GC and MC for Tensor Approach TA (Taubin, 1995)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Initialise values
nrErrorsRcond=0;
% ANGULAR DEFECT
V_AD_new=zeros(nr_vf, 1);  
V_AD_check2=zeros(nr_vf,1);

% Euler's Theorem ET (Watanabe & Belyaev, 2001)
V_GC_ET=zeros(nr_vf,1); 
V_MC_ET=zeros(nr_vf,1);

% Tensor Approach TA (Taubin, 1995)
V_GC_TA=zeros(nr_vf,1);
V_MC_TA=zeros(nr_vf,1);
 
% Errors
V_TT_Err=zeros(nr_vf,3);


%% The main loop; each vertex is processed
VREZ=V;
for i=1:nr_vf % process vertex i
    if V_convexhull (i,1)==0  % the vertex is NOT situated on the convex hull
       varAux=V_1R_NeighNr(i,1);
       mat_v=transpose(V_1R_Neigh(i,1:varAux));
       mat_neigh=V_1R_Neigh_Tri(i,1:varAux);
       [x_vf, y_vf, z_vf] = determineMatrixCoord(V, i, mat_v);  % construct the matrix of coordinates of the neighbours 
       [sum_angles, gaussian_wb, mean_wb, ~, ~, gaussian_t, mean_t, MatErr]=...
           determineValuesETTA(V_weighted_normal, T_area, x_vf, y_vf, z_vf,i, mat_neigh);
       V_AD_new(i,1)=2*pi-sum_angles;
       V_AD_check2(i,1)=V_AD(i,1)-V_AD_new(i,1);
       % Euler's Theorem
       V_GC_ET (i,1)=gaussian_wb;
       V_MC_ET (i,1)=mean_wb;
       % Tensor Approach
       V_GC_TA (i,1)= gaussian_t;
       V_MC_TA (i,1)= mean_t;
       % Errors
       V_TT_Err(i,:)=MatErr(1,:);
    end
end
