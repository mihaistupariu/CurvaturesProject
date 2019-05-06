function [ coord1_vf, coord2_vf, coord3_vf ] = ...
    determineMatrixCoord(vertices, index_vf, matr_coord_vf )

%% Description
% The function generates the coordinates of the 1-ring neigbourhood of
% index_vf RELATIVE TO THIS VERTEX (IT IS THE ORIGIN) 
% Input. 
%   vertices: the matrix of coordinates of the vertices
%   index_vf: the vertex processed
%   matr_coord_vf: the IDs of the neighbours
% Output.
%   coord1_vf, coord2_vf, coord3_vf: column matrices for the x,y,z (1,2,3)
%   coordinates

 
[dim_matr,~]=size (matr_coord_vf); % the size of the matrix, i.e. # of neighbours
coord1_vf=zeros (dim_matr,1); coord2_vf=zeros (dim_matr,1); coord3_vf=zeros (dim_matr,1); 

    for jj=1: dim_matr % take all the neighbours / points
        % find the index of the neighbour
        index_vf_neighh=matr_coord_vf(jj); 
        %translate the coordinates of the point by the ones of the
        %'central' vertex
        coord1_vf(jj,1)=vertices(index_vf_neighh,1)-vertices(index_vf,1); 
        coord2_vf(jj,1)=vertices(index_vf_neighh,2)-vertices(index_vf,2);
        coord3_vf(jj,1)=vertices(index_vf_neighh,3)-vertices(index_vf,3);   
    end
end
 