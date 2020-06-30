function [pc_Output, reg_Grid]= generateFinalPointCloud (pc_Input, cell_Size, x_Min, x_Max, y_Min, y_Max)


%% Description
% Transforms a point cloud into a new point cloud corresponding to a 
% regular grid of a given cell size. In the new point cloud, each cell has
% assigned a single point and conversely. The new points are obtaining by
% averaging in an appropriate manner the original ones.
% Input.
%   pc_Input: the input point cloud
%   cell_Size: the size of the cell
%	x_Min,....: the limits of the point cloud 
% Output.
%   pc_Output: the output point cloud (without (x,y) duplicates)
%   reg_Grid: the elevations in the corresponding regular grid. 
%   Each point of pc_Output corresponds to a single cell of reg_Grid and
%   conversely.

%% Main steps

% Transform in raster with cell size cellSize and get correspondence
[pcIntermediate, ~, reg_Grid,ptsInCell, ~, ~,x_Min_Final,~, ~, y_Max_Final]=...
    linkPointCloudToRaster(pc_Input, cell_Size, x_Min, x_Max, y_Min, y_Max);

% Find additional points needed and fill the gaps
meanVal=mean(pc_Input(:,3));
[additionalPoints, nrAddPoints]=fillGaps(reg_Grid, ptsInCell, x_Min_Final,y_Max_Final, cell_Size,meanVal);
[nrVfAux,~]=size(pcIntermediate);
pcIntermediate(nrVfAux+1:nrVfAux+nrAddPoints,1:3)=additionalPoints(1:nrAddPoints, 1:3);

% Generate again raster / pc
[pcFinal, ~, reg_Grid,~, pcToGrid, ~, ~,~,~,~]=...
    linkPointCloudToRaster(pcIntermediate, cell_Size, x_Min, x_Max, y_Min, y_Max);
 a1=pcFinal;
 a2=pcToGrid;
% Transform the point cloud in the final point cloud
[nrows,ncols]=size(reg_Grid);
v1=computeAveragePC(nrows, ncols, pcToGrid,pcFinal(:,1));
v2=computeAveragePC(nrows, ncols, pcToGrid,pcFinal(:,2));
v3=computeAveragePC(nrows, ncols, pcToGrid,pcFinal(:,3));
pc_Output= [matrixToColumnVector(v1), matrixToColumnVector(v2), matrixToColumnVector(v3)];
%disp(pc_Output);
 
end


        %% FUNCTIONS USED 
        
        %% Fill gaps
        function [pts_added, nr_pts]= fillGaps (matr_in, mat_points_in_cell, xleft, ytop, lung, val_mean)

        %% Description: adds point where needed


        %% Find how many points are needed
        [nr, nc]=size(mat_points_in_cell);
        nr_pts=0;
        for ii=1:nr
            for jj=1:nc
                if mat_points_in_cell(ii,jj)==0
                    nr_pts=nr_pts+1; 
                end
            end
        end

        %% Create the matrix

        pts_added=zeros(nr_pts,3);

        %% Main loop; fill the gaps. ASSUMPTIONS: (i) 'isloated' cells; (ii) not on the border

        index=0;
        for ii=1:nr
            for jj=1:nc
                if mat_points_in_cell(ii,jj)==0
                    index=index+1;
                    x=xleft+(jj-1)*lung+lung/2;
                    y=ytop-(ii)*lung+lung/2;
                    %{
                    z=mat_grid(ii-1,jj-1)+mat_grid(ii-1,jj)+mat_grid(ii-1,jj+1)...
                        +mat_grid(ii,jj-1)+mat_grid(ii,jj+1)+mat_grid(ii+1,jj-1)...
                        +mat_grid(ii+1,jj)+mat_grid(ii+1,jj+1); z=z/8; % elevation is mean of the surrounding cells
                    %}
                    
                    %%%%
                    
                                            %% Clip the "neighbouring" matrices
                                delta=1; % for tests, one could consider larger neighbourhouds

                                % (i) Outside borders
                                    if (ii>1) && (ii<nr) && (jj>1) && (jj<nc)  
                                        matr_loc=matr_in(ii-delta:ii+delta,jj-delta:jj+delta);
                                    end

                                % (ii) The borders without corners
                                    % Upper border, r=1
                                    if ii==1 && jj~=1 && jj~=nc
                                       matr_loc=matr_in(ii:ii+delta,jj-delta:jj+delta);
                                    end
                                    % Lower border, r=nr
                                    if ii==nr && jj~=1 && jj~=nc
                                       matr_loc=matr_in(ii-delta:ii,jj-delta:jj+delta);
                                    end
                                    % Left border, c=1
                                    if jj==1 && ii~=1 && ii~=nr
                                       matr_loc=matr_in(ii-delta:ii+delta,jj:jj+delta);
                                    end
                                    % Right border, c=nc
                                    if jj==nc && ii~=1 && ii~=nr
                                       matr_loc=matr_in(ii-delta:ii+delta,jj-delta:jj);
                                    end     

                                % (iii) The corners
                                    % Up left
                                    if ii==1 && jj==1
                                       matr_loc=matr_in(ii:ii+delta,jj:jj+delta);
                                    end
                                    % Up right
                                    if ii==1 && jj==nc
                                       matr_loc=matr_in(ii:ii+delta,jj-delta:jj);
                                    end
                                    % Down left
                                    if ii==nr && jj==1
                                       matr_loc=matr_in(ii-delta:ii,jj:jj+delta);
                                    end
                                    % Down right
                                    if ii==nr && jj==nc
                                       matr_loc=matr_in(ii-delta:ii,jj-delta:jj);
                                    end 

                            %% If the central point already has a value
                            if matr_in(ii,jj)~=-9999
                                z=matr_in(ii,jj);
                            end

                            %% If the central point is "NODATA" 
                            if matr_in(ii,jj)==-9999
                                % Start analysing the local neighbours
                                sumaLocal=0;
                                indexLocal=0;
                                [randuri,coloane]=size(matr_loc);
                                for rr=1:randuri
                                    for cc=1:coloane
                                        if matr_loc(rr,cc)~=-9999
                                            indexLocal=indexLocal+1;
                                            sumaLocal=sumaLocal+matr_loc(rr,cc);
                                        end
                                    end
                                end
                                if indexLocal~=0
                                    z=sumaLocal/indexLocal;
                                else
                                    z=val_mean;
                                end
                            end
                    %%%%
                    pts_added(index,1)=x;
                    pts_added(index,2)=y;
                    pts_added(index,3)=z;  
                end
            end
        end

        end
        
        %% Compute the average cloud
        
        function [mat_out] = computeAveragePC(nr, nc, pc_poscell, vect_in)

        %% Description
        % Info transfer between PC and grid; the values for PC are transferred to
        % grid via averaging
        % Input: nr, nc (grid charact.). pc_poscell (relation PC <-> grid); vect_in
        % (charact per PC); a COLUMN vector
        % Output: mat_out (charct. per grid cell)
  


        %% Initialization

        mat_out=zeros(nr,nc);
        mat_out_aux=zeros(nr,nc);
        index_aux=zeros(nr,nc);
        [nr_vf, ~]=size(vect_in);


        %% Main loop


        for ii=1:nr_vf
            aa=isnan (vect_in(ii,1));
            if aa==0
                rr=pc_poscell(ii,1); % row
                cc=pc_poscell(ii,2);  %column
                index_aux(rr,cc)=index_aux(rr,cc)+1; % a new element was found
                mat_out_aux(rr,cc)=mat_out_aux(rr,cc)+vect_in(ii,1); % add the value to the cell
            end
        end

        %% Final average

        for ii=1:nr
            for jj=1:nc
               if index_aux(ii,jj)~=0 
                   mat_out(ii,jj)=mat_out_aux(ii,jj)/index_aux(ii,jj); % average
               else 
                   mat_out(ii,jj)=0/0;  % not available
               end 
            end
        end



        end


        %% Transform matrix to column
        function [ vect_out ] = matrixToColumnVector(mat_in)

        %% Description: Transforms a matrix into a vector for statistics
        % Input: the matrix mat_in
        % Output: a column vector containing the elements of the matrix
         


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
