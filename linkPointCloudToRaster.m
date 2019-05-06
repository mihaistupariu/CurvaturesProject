 function [mat_pc_final, mat_pc_exclude, mat_grid , mat_points_in_cell, pc_poscell, grid_pts, x_min_final, x_max_final, y_min_final, y_max_final] =...
            linkPointCloudToRaster(mat_pc, cell_dim, x_min, x_max, y_min, y_max)

%% Description 
%
% Input: mat_pc: coordinates of points (rows: observations, columns: 3 coordinates)
%        cell_dim: size of the cell; 
%        x_min, ... : limits indicated by the user-to be compared with
%        the limits of the point cloud
%        
% Output: mat_pc_final, mat_pc_exclude: new pc (final/excluded)
%         mat_grid: regular grid of elevations (DEM) with resolution cell_dim
%         of a... cell is computed by averaging the values of the points inside the cloud
%         mat_points_in_cell: number of points/cell 
%         pc_poscell: indicates the cell of the grid for each point
%         grid_pts: indicates for each cell the points inside it
%         x_min_final,... : limits of the cropped point cloud
%
% ASSUMPTION: very high point density
 

    %% Initialise the grid
    
    ncols=0;     
    aux=max(x_min, min (mat_pc(:,1))); % compare limit of the point cloud with the indicated one
    x_min_final=floor(aux); % get integer value
    aux=min(x_max, max (mat_pc(:,1))); 
    x_max_final=floor(aux)+1;
    if x_max_final>=x_min_final+1
        ncols=floor((x_max_final-x_min_final)/cell_dim); % number of columns, if it makes sense
    end 
    
    nrows=0;
    aux=max(y_min, min(mat_pc(:,2)));  % compare limit of the point cloud with the indicated one
    y_min_final=floor(aux);  % get integer value
    aux=min(y_max, max(mat_pc(:,2)));
    y_max_final=floor(aux)+1;
    if y_max_final>=y_min_final+1
        nrows=floor((y_max_final-y_min_final)/cell_dim);  % number of rows, if it makes sense
    end

    mat_grid=zeros(nrows, ncols);
    mat_points_in_cell=zeros(nrows, ncols);
    [nr_pt,~]=size(mat_pc);
    pc_poscell=zeros(nr_pt,2);
    mat_pc_final=zeros(nr_pt,3);
    mat_pc_exclude=zeros(nr_pt,3);
    grid_pts=zeros(nrows, ncols);
    index=zeros(nrows,ncols);
    index_include=0;
    index_exclude=0;

    %% PROCESS THE DATA IF IT MAKES SENSE
    
    if nrows~=0 && ncols~=0 
        %% Step 1: run over all points; find for each point the appropriate cell
        for vv=1:nr_pt
            % Compute the difference
            xpt= mat_pc(vv,1);
            ypt= mat_pc(vv,2);
            deltax=floor((xpt-x_min_final)/cell_dim)+1;
            deltay=floor((ypt-y_min_final)/cell_dim);
            if deltax>=1 && deltax<=ncols && deltay>=0 && deltay<=nrows-1 % if the point is inside the area of interest
                mat_grid(nrows-deltay,deltax)=mat_grid(nrows-deltay,deltax)+mat_pc(vv,3);
                mat_points_in_cell(nrows-deltay,deltax)=mat_points_in_cell(nrows-deltay,deltax)+1;
                index_include=index_include+1;
                mat_pc_final(index_include,1:3)=mat_pc(vv,1:3);
            else
                index_exclude=index_exclude+1;
                mat_pc_exclude(index_exclude,1:3)=mat_pc(vv,1:3);
            end
        end
        mat_pc_final=mat_pc_final(1:index_include, 1:3);
        mat_pc_exclude=mat_pc_exclude(1:index_exclude, 1:3);
       

        %% Step 2: run over all grid cells; average the elevation and get the grid values
        for rr=1:nrows
            for cc=1:ncols
                if mat_points_in_cell(rr,cc)>0 
                    mat_grid(rr,cc)=mat_grid(rr,cc)/mat_points_in_cell(rr,cc);
                else
                    mat_grid(rr,cc)=-9999;
                end
            end
        end 
        
        %% Step 3: find the correspondence between the final point cloud and the grid
        % Initializations
        maxNrPts=max(max(mat_points_in_cell));
        grid_pts=zeros(nrows, ncols, maxNrPts);
 
 
        % Main loop
        for vv=1:index_include % take each point from the selected point cloud (mat_pc_include)
            % Find the position in grid 
            xpt= mat_pc_final(vv,1);
            ypt= mat_pc_final(vv,2);
            deltax=floor((xpt-x_min_final)/cell_dim)+1;
            deltay=floor((ypt-y_min_final)/cell_dim);
            row=nrows-deltay;
            column=deltax;
            pc_poscell(vv,1)=row;
            pc_poscell(vv,2)=column;
            % A new point was found in the cell; record its ID
            index(row, column)=index(row, column)+1;
            aux=index(row, column);
            grid_pts(row, column, aux)=vv;   
        end

        %% Check consistency
        for rr=1:nrows
            for cc=1:ncols
              % for each cell the number of points must be the same in the two
              % approaches
              if index(rr,cc)~=mat_points_in_cell(rr,cc)
                  disp('EROARE!!!!')
              end
            end
        end

        
    end


    end