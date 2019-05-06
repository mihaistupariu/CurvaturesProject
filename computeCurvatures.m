function [mat_Gaussian_curvatures, mat_Mean_curvatures, vert_CH] = computeCurvatures(V, alpha_vector, cellSize, regGrid)

%% Description
% The function provides the values of Gaussian and mean curvatures for a
% point cloud V corresponding to a regular grid regGrid (it is assumed that
% each point of the point cloud corresponds to a cell of the grid)
% Input.
%   V: the input (original) point cloud
%   alpha_vector: the vector of alpha-values (Integral Approach)
%   cellSize: the size of the cell for which the grid is generated
%   regGrid: the regular grid
% Output.
%   mat_Gaussian_curvatures: matrix of Gaussian curvatures (each column
%   corresponds to one of the methods)
%   mat_Mean_curvatures: matrix of mean curvatures (each column
%   corresponds to one of the methods)
%   vert_CH: vertices situated on the border of the convex hull


%% Main steps - the computations are performed in the scripts
% get topological information
script_0
vert_CH=V_convexhull;

% compute descriptors related to the triangles
script_1

% compute descriptors related to the edges
script_2

% compute the curvatures - Gauss-Bonnet schemes
script_3_GB

% compute the curvatures - Euler Theorem and Tensor Approach
script_4_ET_TA

% compute the curvatures - Jet Fitting
script_5_JF

% compute the curvatures - Normal Cycle
script_6_NC

% compute the curvatures - Integral Approach
script_7_IA

% compute the curvatures - Regular Grid
[A_GC_RG, A_MC_RG]=determineDEMCurvatures(regGrid, cellSize);
V_GC_RG=matrixToColumnVector(A_GC_RG);
V_MC_RG=matrixToColumnVector(A_MC_RG);


%% Gather all vectors of curvatures into the corresponding matrices
% Gaussian curvatures
mat_Gaussian_curvatures=[V_GC_1 V_GC_2 V_GC_ET V_GC_TA V_GC_JF V_GC_NC_1R V_GC_NC_2R V_GC_IA V_GC_RG];
% Mean curvatures
mat_Mean_curvatures=[V_MC_1 V_MC_2 V_MC_ET V_MC_TA V_MC_JF V_MC_NC_1R V_MC_NC_2R V_MC_IA V_MC_RG];

end

%% FUNCTIONS USED
      %% Transform matrix to column
        function [ vect_out ] = matrixToColumnVector(mat_in)

        %% Description: 
        % Transforms a matrix into a vector for statistics
        % Input.
        %   mat_in: the input matrix
        % Output.
        %   vect_out: a column vector containing the elements of the matrix

        
        %% Initialization 
        [nr,nc]=size(mat_in);
        vect_out=zeros(nr*nc,1);

        %% Main loop: transfer elements of the matrix to a vector
        index=0;
        for ii=1:nr
            for jj=1:nc
             index=index+1;
             vect_out(index)=mat_in(ii,jj);
            end
        end

        end
