function [vect_output1, vect_output2] = computeCorrelationVector(vect_ref, mat_in, vect_border)


%% Description
% This function determines the correlations between a given vector of curvatures and
% other curvatures
% Input.
%   vect_ref: the column corresponding to the reference curvature
%   mat_in: the matrix of curvatures - each column represents the values of
%   a curvature computed with a certain method
%   vect_border: column vector indicating with an 1 the vertices on the
%   convex hull; these vertices are excluded from the statistics
% Output.
%   vect_output1: the matrix of correlations between the values of different 
%           curvatures,as provided by the methods considered, at different 
%           levels of resolution. The rows correspond to the
%           levels of resolution. The columns correspond to the methods compared
%           with the reference method. 
%   vect_outputw: the matrix of absolute errors the values of different 
%           curvatures,as provided by the methods considered, at different 
%           levels of resolution. The rows correspond to the
%           levels of resolution. The columns correspond to the methods compared
%           with the reference method. 
 

%% Exclude outliers
[nr_vf,nMethods]=size(mat_in);
vect_exclude=zeros(nr_vf,1);
vect_output1=zeros(1,nMethods);
vect_output2=zeros(1,nMethods);

%% Determine the vertices to be excluded, according to suitable criteria

for ii=1:nr_vf
       % Criterion 1: -9999 (NODATA) for the grid
       if mat_in(ii,nMethods)==-9999
           vect_exclude(ii,1)=1;
       end
       % Criterion 2: on the border
       if vect_border(ii,1)==1
           vect_exclude(ii,1)=1;
       end
end


%% Generate matrices for computing correlation

nrows=nr_vf-sum(vect_exclude);
mat_aux=zeros(nrows, nMethods);
vect_aux=zeros(nrows,1);
index=0;

for ii=1:nr_vf
    if vect_exclude(ii,1)==0
        index=index+1;
        mat_aux(index,1:nMethods)=mat_in(ii,1:nMethods);
        vect_aux(index,1)=vect_ref(ii,1);
    end
end

%% Compute the vector of correlations / absolute errors

for mm=1:nMethods
    auxCorr=corrcoef(mat_aux(:,mm), vect_aux);
    aux=auxCorr(1,2);
    vect_output1(1,mm)=aux;
    vect_output2(1,mm)=averageError(mat_aux(:,mm)-vect_aux);
end

end


%% FUNCTION USED 

    function [ avgerr ] = averageError( vect_column)

    [nr, ~]=size(vect_column);
    avgerr=0;

    for ii=1:nr
        avgerr=avgerr+abs(vect_column(ii,1));
    end

    avgerr=avgerr/nr;
    end
