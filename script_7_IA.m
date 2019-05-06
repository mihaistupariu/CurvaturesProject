%% SCRIPT_7_IA: 

%% Description
% Implementation of the Integral Approach (Pottmann et al., 2007) 
% Computes the 3x3 covariance matrix, then extracts the eigenvalues and
% provides the principal & Gaussian/mean curvatures
% The script uses the alpha-vector read from the GUI

%% Initialization
refRadius=zeros(nr_vf,1); 
[~,nr_alphas]=size(alpha_vector);
V_GC_IA=zeros(nr_vf,nr_alphas);
V_MC_IA=zeros(nr_vf,nr_alphas);

%% Take each coefficient radius from the vector of possible coefficients 

for aa=1:nr_alphas
    coeffRadius=alpha_vector(aa);

        %% Compute the reference radius; it depends on the vertices
        for vv=1:nr_vf %take all vertices
            PctCen=V(vv,1:3); % the vertex vv is central
            matAux=setdiff(V_2R_Neigh(vv,:),0); % take the vertices is the 2-RN and remove 0s
            [~,ncAux]=size(matAux);
            MatDist=zeros(ncAux,1); %matrix for distances
            for ii=1:ncAux 
                ww=V_2R_Neigh(vv,ii); % extract the vertex index on position ii
                PctVar=V(ww,1:3); % get the coordinates
                MatDist(ii,1)=norm(PctVar-PctCen); 
            end
            minAux=min(MatDist);
            refRadius(vv,1)=minAux;
        end

        %% Main loop: process all vertices
        for vv=1:nr_vf
            %% The radius and other initializations
            crtRadius=refRadius(vv,1)*coeffRadius;
            matCova=zeros(3,3); % For each vertex this matrix is initialized as zero
            vfCtrl=V(vv,1:3);
            index=0;
            ariiTri=[];
            baricentreTri=[];
            normalaZonei=[0 0 0];

            %% Process the triangles in the 1-RingN
            set_tri=setdiff(V_1R_Neigh_Tri(vv,:),0);
            [~,numar_tri]=size(set_tri);
            for nn=1:numar_tri %take each triangle
               tri_crt=set_tri(1,nn); % transfer index in the vector to triangle index
               ind1=T(tri_crt,1); ind2=T(tri_crt,2); ind3=T(tri_crt,3);  % indices of the vertices for the triangle tri_crt
               vf1=V(ind1, 1:3); vf2=V(ind2, 1:3); vf3=V(ind3, 1:3); % coordinates of the vertices 
               % Main computation: add contribution to the covariance matrix for
               % the vertex vv
               matCova=matCova+determineCovarianceMatrixTriangleGeneral (vfCtrl, crtRadius, vf1, vf2, vf3);
               % Compute the main measures of the triangle
               [areaT, barT] = determineMeasuresTriangle(vf1, vf2, vf3);  % area and barycenter
               index=index+1; % a new triangle was found
               ariiTri(index,1)=areaT; % add the area to the vector of areas
               baricentreTri(index,1:3)=barT; % add the barycenter to the vector of barycenters
               normalaZonei=normalaZonei+T_normal(tri_crt,1:3); % add the normal contribution
            end

            %% Process the triangles in the 2-RingN
            set_tri=setdiff(V_2R_Neigh_Tri(vv,:),0);
            [~,numar_tri]=size(set_tri);
            for nn=1:numar_tri %take each triangle
               tri_crt=set_tri(1,nn); % transfer index in the vector to triangle index
               ind1=T(tri_crt,1); ind2=T(tri_crt,2); ind3=T(tri_crt,3);  % indices of the vertices for the triangle tri_crt
               vf1=V(ind1, 1:3); vf2=V(ind2, 1:3); vf3=V(ind3, 1:3); % coordinates of the vertices
               % Main computation: add contribution to the covariance matrix for
               % the vertex vv
               matCova=matCova+determineCovarianceMatrixTriangleGeneral (vfCtrl, crtRadius, vf1, vf2, vf3);   
               % Compute the main measures of the triangle
               [areaT, barT] = determineMeasuresTriangle(vf1, vf2, vf3);  % area and barycenter
               index=index+1; % a new triangle was found
               ariiTri(index,1)=areaT; % add the area to the vector of areas
               baricentreTri(index,1:3)=barT; % add the barycenter to the vector of barycenters
               normalaZonei=normalaZonei+T_normal(tri_crt,1:3); % add the normal contribution
            end

            %% Determine the eigenvalues of the covariance matrix
            eigenMatCov=eig(matCova);
            eigVal1=eigenMatCov(1,1); eigVal2=eigenMatCov(2,1); eigVal3=eigenMatCov(3,1); 

            %% Compute the Gaussian and the mean curvature
            % Gaussian curvature
            GaussCurv=36/(5*crtRadius^2)-24/(5*pi*crtRadius^6)*(3*eigVal1+3*eigVal2+2*eigVal3);
            V_GC_IA(vv,aa)=GaussCurv;
            % Absolute value of mean curvature
            absMeanCurv=sqrt(abs(16*eigVal3/(pi*crtRadius^6)+2*GaussCurv/3));
            % Sign of mean curvature
            areaTot=sum(ariiTri); % total area of the 2-R neigborhood
            baricentru=[0 0 0];  % initialize computation for barycenter
            for indexAux=1:index
                baricentru=baricentru+ariiTri(indexAux,1)*baricentreTri(indexAux,1:3);
            end
            baricentru=baricentru/areaTot; % coordinates of the barycenter
            if norm(normalaZonei)~=0
                normala=normalaZonei/norm(normalaZonei); % normalize the sum over the contributors
            else
                normala=[0 0 1];
            end
            normalaInt=-normala; % The interior normal (theoretically downwards)
            semnMeanCurv=sign(dot(baricentru-V(vv,1:3),normalaInt));  % If the vector from the center V(vv,:) to the barycenter is 'interior', the sign is +, otherwise it is -
            V_MC_IA(vv,aa)=semnMeanCurv*absMeanCurv;

        end 
end
