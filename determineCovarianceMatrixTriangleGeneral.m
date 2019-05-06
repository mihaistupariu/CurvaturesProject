function [mat_cova] = determineCovarianceMatrixTriangleGeneral (vfCentral, raza, vf_A, vf_B, vf_C)

%% Description 
% Computes the covariance matrix for the intersection between the spehere
% of center vfCentral and radius raza
% the triangle determined by the vertices vfA, vfB, vfC
%  

%% Initialization
% Translate such that the origin is at vfCentral

vfA=vf_A-vfCentral; normA=norm(vfA);
vfB=vf_B-vfCentral; normB=norm(vfB);
vfC=vf_C-vfCentral; normC=norm(vfC);

%% Preparations
% If the vectors are line vectors, transform them into column vectors
[nr,nc]=size(vfA);
if nr==1 && nc==3 
    vfA=transpose(vfA);
end

[nr,nc]=size(vfB);
if nr==1 && nc==3 
    vfB=transpose(vfB);
end

[nr,nc]=size(vfC);
if nr==1 && nc==3 
    vfC=transpose(vfC);
end

%% All points inside the ball
if normA<=raza && normB<=raza && normC<=raza 
    mat_cova=determineCovarianceMatrixTriangle(vfA, vfB, vfC);
end


%% Two points inside the ball
if normA<=raza && normB<=raza && normC>raza  % A and B inside; C outside
    pA=determineIntersectionSphereSegment(raza, vfA, vfC );
    pB=determineIntersectionSphereSegment(raza, vfB, vfC );
    mat_cova=determineCovarianceMatrixTriangle(vfA, pA, pB)+determineCovarianceMatrixTriangle(vfB, pA, pB);
end

if normB<=raza && normC<=raza && normA>raza  % B and C inside; A outside
    pB=determineIntersectionSphereSegment(raza, vfB, vfA );
    pC=determineIntersectionSphereSegment(raza, vfC, vfA );
    mat_cova=determineCovarianceMatrixTriangle(vfB, pB, pC)+determineCovarianceMatrixTriangle(vfC, pB, pC);
end

if normC<=raza && normA<=raza && normB>raza  % C and A inside; B outside
    pC=determineIntersectionSphereSegment(raza, vfC, vfB );
    pA=determineIntersectionSphereSegment(raza, vfA, vfB );
    mat_cova=determineCovarianceMatrixTriangle(vfC, pC, pA)+determineCovarianceMatrixTriangle(vfA, pC, pA);
end

%% One point inside the ball
if normA<=raza && normB>raza && normC>raza  % A inside; B, C outside
    pB=determineIntersectionSphereSegment(raza, vfA, vfB );
    pC=determineIntersectionSphereSegment(raza, vfA, vfC );
    mat_cova=determineCovarianceMatrixTriangle(vfA, pB, pC);
end

if normB<=raza && normC>raza && normA>raza  % B inside; C, A outside
    pC=determineIntersectionSphereSegment(raza, vfB, vfC );
    pA=determineIntersectionSphereSegment(raza, vfB, vfA );
    mat_cova=determineCovarianceMatrixTriangle(vfB, pC, pA);
end

if normC<=raza && normA>raza && normB>raza  % C inside; A, B outside
    pA=determineIntersectionSphereSegment(raza, vfC, vfA );
    pB=determineIntersectionSphereSegment(raza, vfC, vfB );
    mat_cova=determineCovarianceMatrixTriangle(vfC, pA, pB);
end
%% All points outside the ball
if normA>raza && normB>raza && normC>raza 
    mat_cova=zeros(3,3);
end



end

    %% FUNCTIONS USED

    function [mat_cova ] = determineCovarianceMatrixTriangle( A, B, C )

    %% Description
    % Implementation of the covariance matrix - Pottmann et al. 2007

    %% Preparations
    % If the vectors are line vectors, transform them into column vectors
    [nr,nc]=size(A);
    if nr==1 && nc==3 
        A=transpose(A);
    end

    [nr,nc]=size(B);
    if nr==1 && nc==3 
        B=transpose(B);
    end

    [nr,nc]=size(C);
    if nr==1 && nc==3 
        C=transpose(C);
    end

    % Coordinates of the points
    xA=A(1,1); yA=A(2,1); zA=A(3,1);
    xB=B(1,1); yB=B(2,1); zB=B(3,1);
    xC=C(1,1); yC=C(2,1); zC=C(3,1);

    x=[xA, xB, xC];
    y=[yA, yB, yC];
    z=[zA, zB, zC];



    aux(1,1)=(dot(x,x)+sum(x)*sum(x))/24;
    aux(1,2)=(dot(x,y)+sum(x)*sum(y))/24; 
    aux(1,3)=(dot(x,z)+sum(x)*sum(z))/24; 
    aux(2,1)=(dot(y,x)+sum(y)*sum(x))/24;
    aux(2,2)=(dot(y,y)+sum(y)*sum(y))/24;
    aux(2,3)=(dot(y,z)+sum(y)*sum(z))/24;
    aux(3,1)=(dot(z,x)+sum(z)*sum(x))/24;
    aux(3,2)=(dot(z,y)+sum(z)*sum(y))/24; 
    aux(3,3)=(dot(z,z)+sum(z)*sum(z))/24; 

    %disp(aux)
    %disp(rank(aux))

    [area, bar] = determineMeasuresTriangle(A, B, C);
    mat_cova=aux-area*bar*transpose(bar);

    end



