function [gauss_curv, mean_curv] = determineCurvaturesJetFitting(A10, A01, A20, A11, A02)

%% Description 
% Computes the curvature for an explicit polynomial surface
% z=A00+A10x+A01y+1/2(A20x^2+2A11xy+A20y^2)


firstForm=[1+A10^2, A10*A01; A10*A01, 1+A01^2];
normainv=1/(sqrt(1+A10^2+A01^2));
secondForm=[normainv*A20, normainv*A11; normainv*A11, normainv*A02];
opWeingarten=firstForm\secondForm;
eigenOpW=eig(opWeingarten);
eigVal1=eigenOpW(1,1);
eigVal2=eigenOpW(2,1);
gauss_curv=eigVal1*eigVal2;
mean_curv=(eigVal1+eigVal2)/2;

end

