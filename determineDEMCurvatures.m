function [mat_GC, mat_MC] = determineDEMCurvatures(mat_elev, lung)

%% Description
% Computes the curvature for a regular grid == DEM 
% Input.
%   mat_elev: matrix of elevations
%   lung: CellSize
% Output.
%   mat_GC, mat_MC: matrices of Gaussian/mean Curvatures
 

%% Initialization
[nrows, ncols]=size(mat_elev);
mat_GC=zeros(nrows, ncols);
mat_MC=zeros(nrows, ncols);

% Put NODATA on the borders
for ir=1:nrows
    mat_GC(ir,1)=-9999; mat_GC(ir,ncols)=-9999;
    mat_MC(ir,1)=-9999; mat_MC(ir,ncols)=-9999;
end
for ic=1:ncols
    mat_GC(1,ic)=-9999; mat_GC(nrows,ic)=-9999;
    mat_MC(1,ic)=-9999; mat_MC(nrows,ic)=-9999;
end

%% Main loop: compute the value for each 'valid' cell

for ir=2:nrows-1 % for loop for rows; from 2:nrows-1; avoid the border
    for ic=2:ncols-1 % for loop for columns; same rule
          % all computations make sense if ~=NODATA  and if the neighbors must be also ~=NODATA
            if (mat_elev(ir,ic)~=-9999) && (mat_elev(ir-1,ic-1)~=-9999) && (mat_elev(ir-1,ic)~=-9999) && (mat_elev(ir-1,ic+1)~=-9999) && (mat_elev(ir,ic-1)~=-9999) && (mat_elev(ir,ic+1)~=-9999) && (mat_elev(ir+1,ic-1)~=-9999) && (mat_elev(ir+1,ic)~=-9999) && (mat_elev(ir+1,ic+1)~=-9999)
                % the cell is "valid"
                % elevations of the cells in the 3x3 moving window - 
                %  notation from Hengl
                z1= mat_elev(ir-1,ic-1);
                z2= mat_elev(ir-1,ic);
                z3= mat_elev(ir-1,ic+1);
                z4= mat_elev(ir,ic-1);
                z5= mat_elev(ir,ic);
                z6= mat_elev(ir,ic+1);
                z7= mat_elev(ir+1,ic-1);
                z8= mat_elev(ir+1,ic);
                z9= mat_elev(ir+1,ic+1);
                
                
                % compute the 'partial derivatives' of the height map zeta
                % (method= quadric fitting, see Y-E, Pennock, Hengl, Shary)
                dfdx=(z3+z6+z9-z1-z4-z7)/(6*lung);
                dfdy=(z1+z2+z3-z7-z8-z9)/(6*lung);
                d2fdx2=(z3+z6+z9+z1+z4+z7-2*z2-2*z5-2*z8)/(3*lung*lung);
                d2fdy2=(z1+z2+z3+z7+z8+z9-2*z4-2*z5-2*z6)/(3*lung*lung);
                d2fdxdy=(z3+z7-z1-z9)/(4*lung*lung);

                % use the equivalent formulas (as in Shary)
                norma=1+dfdx*dfdx+dfdy*dfdy;
                %%%
                numarator_kg=d2fdx2*d2fdy2-d2fdxdy*d2fdxdy;
                numitor_kg=norma.^2;
                mat_GC(ir,ic)=numarator_kg/numitor_kg;
                %%%
                numarator_h=(1+dfdx^2)*d2fdy2+(1+dfdy^2)*d2fdx2-2*dfdx*dfdy*d2fdxdy;
                numitor_h=2*(sqrt(norma))^3;
                mat_MC(ir,ic)=numarator_h/numitor_h;  
                
            else
                mat_GC(ir,ic)=-9999; 
                mat_MC(ir,ic)=-9999;
            end
    end
    
end


end

