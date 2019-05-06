function [area, bar] = determineMeasuresTriangle(A, B, C)

%% Description 
% Computes basic measures (area, barycenter) of a triangle ABC
% Input.
%	A, B, C: the vertices
% Output.
%   area, bar: the area and the barycenter of the triangle
 


aux=cross (B-A, C-A);
area=0.5*norm(aux); 
bar=(A+B+C)/3;



end

