function [A00, A10, A01, A20, A11, A02] = determineCoefficientsFP(coordFittedBasis)

%% Description
% Determine coefficients of the fitted polynomial (cf. Cazals & Pouget,
% 2008)

%% Initializations
[~, nrPoints]=size(coordFittedBasis); % this is a column vector

%% Prepare the matrices for the fitting
matZ=transpose(coordFittedBasis(3,1:nrPoints));

matX=transpose(coordFittedBasis(1,1:nrPoints));
matY=transpose(coordFittedBasis(2,1:nrPoints));

% Construct the matrix 
matM(:,1)=ones(nrPoints,1);
matM(:,2)=matX;
matM(:,3)=matY;

for ii=1:nrPoints
   matM(ii,4)=matX(ii,1)*matX(ii,1)/2;
   matM(ii,5)=matX(ii,1)*matY(ii,1);
   matM(ii,6)=matY(ii,1)*matY(ii,1)/2;
end
% Mean value
valh=mean(makeVector(matM(1:nrPoints,2:3)));
if (valh==0)
    valh=0.0000001;
end
% Modify the matrices and adjust by mean value
matD=diag([1 valh valh valh^2 valh^2 valh^2]);
matMprim=matM/matD;
% Singular value decomposition for the matrix matMprim
[matU, matS, matV]=svd(matMprim);
% For the diagonal part, determine inverse
matS_new=zeros(6,nrPoints);
rankS=rank(matS);
for ii=1:rankS
    if matS(ii,ii)~=0
        matS_new(ii,ii)=1/matS(ii,ii);
    end
end
% Find solution for the new system
matY=matV*matS_new*transpose(matU)*matZ;
matA=matD\matY;
% Extract coefficients
A00=matA(1,1);
A10=matA(2,1);
A01=matA(3,1);
A20=matA(4,1);
A11=matA(5,1);
A02=matA(6,1);



end

%% FUNCTIONS USED
function [ vect_out ] = makeVector( matr_in )

[nr, nc]=size(matr_in);
vect_out=zeros(1,nr*nc);
ima=0;
for ii=1:nr
    for jj=1:nc 
        ima=ima+1;
        vect_out(ima)=matr_in(ii,jj);
    end
end
end


