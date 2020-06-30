fid=fopen('0_function.txt');
%
line1=fgetl(fid);
func1=str2func(line1);
f=func1(X,Y);
%
line2=fgetl(fid);
func2=str2func(line2);
dfdx=func2(X,Y);
%
line3=fgetl(fid);
func3=str2func(line3);
dfdy=func3(X,Y);
%
line4=fgetl(fid);
func4=str2func(line4);
d2fdx2=func4(X,Y);
%
line5=fgetl(fid);
func5=str2func(line5);
d2fdxdy=func5(X,Y);
%
line6=fgetl(fid);
func6=str2func(line6);
d2fdy2=func6(X,Y);
%
%% Description
% Implements the computation of Gaussian and mean curvatures, when partial
% derivatives are known and the surface is a height function, i.e. 
% (x,y)|--> (x, y, f(x,y)) 
% Input: Partial derivatives: 1st order dfdx, dfdy; 2nd order: d2fdx2,  d2fdxdy, d2fdy2
% Output Gaussian curvature k and mean curvature h
% 15.01.2015, Mihai Sorin Stupariu


%% Auxiliary: the norm 
norma=1+dfdx.*dfdx+dfdy.*dfdy;

%% Gaussian curvatures

numarator_kg=d2fdx2.*d2fdy2-d2fdxdy.*d2fdxdy;
numitor_kg=norma.^2;
GCC=numarator_kg./numitor_kg;

%% Mean curvature

numarator_h=(1+dfdx.^2).*d2fdy2+(1+dfdy.^2).*d2fdx2-2.*dfdx.*dfdy.*d2fdxdy;
numitor_h=2.*(sqrt(norma)).^3;
MCC=numarator_h./numitor_h;


 