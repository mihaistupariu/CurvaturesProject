function [coordFittedBasis] = determineCoordFittingBasis(P0Row,coordCanBasisRows)

%% Description
% Determine coefficients of the fitting basis (cf. Cazals & Pouget,
% 2008)


%% Initialization
% IT IS ASSUMED THAT THE POINTS ARE IN ROW FORM: TRANSPOSE TO GET COLUMN
% FORM
coordCanBasis=transpose(coordCanBasisRows); % transform the matrix for having the points on the columns
P0=transpose(P0Row);
[~,nrPoints]=size(coordCanBasis);  % number of points
coordCanBasisMinusCG=zeros(3,nrPoints); % initialize the matrix for which the covariance Matrix is determined
coordCanBasisMinusP0=zeros(3,nrPoints); % initialize the matrix after translation
% The barycenter CG
coordCG=[sum(coordCanBasis(1,:))/nrPoints;  sum(coordCanBasis(2,:))/nrPoints; sum(coordCanBasis(3,:))/nrPoints];

%% Compute the 3x3 covariance matrix
% Substract the barycenter CG
for ii=1:nrPoints
    coordCanBasisMinusCG(1:3,ii)=coordCanBasis(1:3,ii)-coordCG; % substract the CG; sum of columns is 0
end
covarianceMatrix=coordCanBasisMinusCG*transpose(coordCanBasisMinusCG);

%% Determine the eigenvalues and the eigenvectors
eigenValues=eig(covarianceMatrix);
[eigenVectors,~]=eig(covarianceMatrix);
% Ordering the eigenvalues / eigenvectors
[minAbsEig,~]=find (abs(eigenValues)==min(abs(eigenValues))); % find the index with eig of smallest amplitude
indiceMinAbsEig=minAbsEig(1,1);
if indiceMinAbsEig~=3 % if not on the last position
    eig_proviz=eigenValues(indiceMinAbsEig);
    eigenValues(indiceMinAbsEig,1)=eigenValues(3,1);
    eigenValues(3,1)=eig_proviz;
    eigenVect_proviz=eigenVectors(:,indiceMinAbsEig);
    eigenVectors(:,indiceMinAbsEig)=eigenVectors(:,3);
    eigenVectors(:,3)=eigenVect_proviz; 
end
% Change the sign of the determinant
if det(eigenVectors)<0
    eigenVectors(:,3)=-eigenVectors(:,3);
end
%% Change of Frame
changeBasis=transpose(eigenVectors); % the matrix that transforms the COLUMN vectors of the eigenbasis in the canonical basis e1 e2 e3

%% New coordinates, in the fitted basis
% Translation such that the central point P0 becomes origin
for ii=1:nrPoints
    coordCanBasisMinusP0(1:3,ii)=coordCanBasis(1:3,ii)-P0; % substract the CG; sum of columns is 0
end
coordFittedBasis=changeBasis*coordCanBasisMinusP0;
