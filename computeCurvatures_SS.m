function [mat_Gaussian_curvatures, mat_Mean_curvatures, vert_CH] = computeCurvatures_SS(V, alpha_vector, cellSize, file_PD )

%% Description
% The function provides the values of Gaussian and mean curvatures for a
% point cloud V corresponding to a smooth parametrized surface (x,y)|---> f(x,y).
% The partial derivatives of the components are read from the file
% 'file_PD'
% Input.
%   V: the input (original) point cloud
%   alpha_vector: the vector of alpha-values (Integral Approach)
%   cellSize: the size of the cell for which the grid is generated
%   file_PD: text file with the partial derivatives of f
% Output.
%   mat_Gaussian_curvatures: matrix of Gaussian curvatures (each column
%   corresponds to one of the methods)
%   mat_Mean_curvatures: matrix of mean curvatures (each column
%   corresponds to one of the methods)
%   vert_CH: vertices situated on the border of the convex hull


%% Main steps 
% The computations are performed in the scripts. The variables (column
% vectors) corresponding to various methods are generated in the scripts

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

% compute the curvatures - Smooth Surface
[V_GC_SS, V_MC_SS]=computeCurvaturesSmoothFunction(V(:,1), V(:,2), file_PD);


%% Gather all vectors of curvatures into the corresponding matrices
% Gaussian curvatures
mat_Gaussian_curvatures=[V_GC_1 V_GC_2 V_GC_ET V_GC_TA V_GC_JF V_GC_NC_1R V_GC_NC_2R V_GC_IA V_GC_SS];
% Mean curvatures
mat_Mean_curvatures=[V_MC_1 V_MC_2 V_MC_ET V_MC_TA V_MC_JF V_MC_NC_1R V_MC_NC_2R V_MC_IA V_MC_SS];

 

end

%% FUNCTIONS USED
        
        %% Computation of Gaussian and mean curvatures
        function [vect_gauss, vect_mean] = computeCurvaturesSmoothFunction(X, Y, file_name)
        %% Description
        % Implements the computation of Gaussian and mean curvatures, when partial
        % derivatives are known and the surface is a height function, i.e. 
        % (x,y)|--> (x, y, f(x,y)) 
        % Input.
        %	X, Y: coordinate vectors
        %   file_name: file with the partial derivatives
        % Output.
        %   vect_gauss, vect_mean: vectors of Gaussian and of mean curvatures
 
        %% Read from file and compute partial derivatives 
        frewind(file_name);
        %
        line1=fgetl(file_name);
        func1=str2func(line1);
        f=func1(X,Y);
        %
        line2=fgetl(file_name);
        func2=str2func(line2);
        dfdx=func2(X,Y);
        %
        line3=fgetl(file_name);
        func3=str2func(line3);
        dfdy=func3(X,Y);
        %
        line4=fgetl(file_name);
        func4=str2func(line4);
        d2fdx2=func4(X,Y);
        %
        line5=fgetl(file_name);
        func5=str2func(line5);
        d2fdxdy=func5(X,Y);
        %
        line6=fgetl(file_name);
        func6=str2func(line6);
        d2fdy2=func6(X,Y);
        
        %% Auxiliary: the norm 
        norma=1+dfdx.*dfdx+dfdy.*dfdy;

        %% Gaussian curvatures
        numarator_kg=d2fdx2.*d2fdy2-d2fdxdy.*d2fdxdy;
        numitor_kg=norma.^2;
        vect_gauss=numarator_kg./numitor_kg;

        %% Mean curvature

        numarator_h=(1+dfdx.^2).*d2fdy2+(1+dfdy.^2).*d2fdx2-2.*dfdx.*dfdy.*d2fdxdy;
        numitor_h=2.*(sqrt(norma)).^3;
        vect_mean=numarator_h./numitor_h;

        end
