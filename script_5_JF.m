%% SCRIPT_5_JF: 

%% Description
% Implementation of the Jet Fitting Approach (Cazals & Pouget, 2008) 

%% Initialization
 
V_GC_JF=zeros(nr_vf,1);
V_MC_JF=zeros(nr_vf,1);

%% Main loop: process all vertices

for vv=1:nr_vf
    % Take all points in the 1-RN and 2-RN
    set_1=setdiff(V_1R_Neigh(vv,:),0);
    set_2=setdiff(V_2R_Neigh(vv,:),0);
    set_VeciniPC=union(set_1,set_2);
    % Number of points 
    [~,nr_VeciniPC]=size(set_VeciniPC);
    coordinatesCanBasisRows=zeros(nr_VeciniPC,3);
    % Coordinates of the central point
    PCen=V(vv,1:3);
    % Extract coordinates of the neighbours
    for ii=1:nr_VeciniPC % take each neighbour
        ww=set_VeciniPC(1,ii);  % transfer to index
        coordinatesCanBasisRows(ii,1:3)=V(ww,1:3); % the point with index ii is on position ww in the original PC V
    end
    % Determine coordinates of the fitted basis (PCA+change of basis)
    [coordinatesFittedBasis] = determineCoordFittingBasis(PCen,coordinatesCanBasisRows);
    % Compute coefficients of the fitting polynomial surface (SVD)
    [a00, a10, a01, a20, a11, a02] =determineCoefficientsFP(coordinatesFittedBasis);
    % Compute curvatures for the fitted surface
    [GCurv, MCurv] = determineCurvaturesJetFitting(a10, a01, a20, a11, a02);
    V_GC_JF(vv,1)=GCurv;
    V_MC_JF(vv,1)=-MCurv;
end
