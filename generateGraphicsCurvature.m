function [  ] = generateGraphicsCurvature (mat_elev, vect_col, az, elev)

%% Description
% Provides the graphical representation for a single curvature. The
% elevations are provided by the regular grid matrix, while the colors are
% provided by the associated curvatures.
% Input.
%   mat_elev: elevations of the regular grid
%   vect_col: vector of curvatures, provides the colors for the cells of
%   the grid
%   az, elev: azimuth, elevation - parameters for 3D view
% Output. The grahical representation

hold off

%% Initializations
[nr, nc]=size(mat_elev);
mat_col=zeros(nr,nc);
 
%% Transfer the color vector to a matrix
mat_aux=ColumnVectorToMatrix(vect_col, nr, nc);
 
%% Replace with zeros on the border
for ii=1:nr
    mat_aux(ii,1)=0; mat_aux(ii,nc)=0;
end
for jj=1:nc
    mat_aux(1,jj)=0; mat_aux(nr,jj)=0;
end

%% Replace 'NODATA' with 0 
for ii=1:nr
    for jj=1:nc
        if mat_aux(ii,jj)==-9999
            mat_col(ii,jj)=0;
        else
            mat_col(ii,jj)=mat_aux (ii,jj);
        end
    end
end
 

%% 3D View
surf (mat_elev,mat_col);
view(az, elev);  
 
end

%% FUNCTION USED
    function [ mat_out] = ColumnVectorToMatrix(vect_in, nr, nc)

    %% Description: Transforms a matrix into a vector for statistics
    % Input: the column vector vect_in
    % Output: a row vector containing the elements of the matrix
 


    %% Initialization 
    [nr_vect,~]=size(vect_in);
    mat_out=zeros(nr,nc);

    %% Main loop: transfer elements of the matrix to a vector

    if nr_vect~=nr*nc 
        disp('Eroare!!!')
    end
    index=0;
    for ii=1:nr
        for jj=1:nc
         index=index+1;
         mat_out(ii,jj)=vect_in(index);
        end
    end

    end

