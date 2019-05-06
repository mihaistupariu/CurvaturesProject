function [pc_Output] = generateAuxPointCloud(pc_Input, x_Min, x_Max, y_Min, y_Max)

%% Description
% Generates a point cloud without (x,y) duplicates
% Input.
%   pc_Input: the input point cloud
%	x_Min,....: the limits of the point cloud 
% Output.
%    pc_Output: the output point cloud (without (x,y) duplicates)
 
%% Main steps
% Transform in raster with cell size 1 and get the correspondence between
% the point cloud and the raster for handling the points
[pcIntermediate, ~, ~,ptsInCell, pcToGrid, gridToPC, ~,~,~,~]=...
    linkPointCloudToRaster(pc_Input, 1, x_Min, x_Max, y_Min, y_Max);
               
% Remove the duplicates by using the regular grid as an auxiliary tool
[pc_Output, index, ~]=removeDuplicates(pcIntermediate, pcToGrid, gridToPC, ptsInCell);
 
end



%% FUNCTIONS USED 
        % removeDuplicates
        function [pc_out, index, index_2] = removeDuplicates(pc_in, pc_poscell, grid_pts, mat_points_in_cell)

        %% Description
        % searches xy_duplicates (rows having same elements on 1st and 2nd positions
        % and counts them
        % INPUT: PC INFO (corresponding cells); GRID INFO (pc_poscell -- IDs of
        %        points, grid_pts, mat_points_in_cell -- nr. of points per cell)
        % OUTPUT: PC without duplicates (pc_out); index (number of duplicates)
         %% Initialization 

        [nr_vf,~]=size(pc_in);
        mat_control=zeros(nr_vf,1);
        index=0;
        pc_in_xy=pc_in(1:nr_vf, 1:2);

        %% Main loop take all rows and search duplicates

        for ii=1:nr_vf 
            if mat_control (ii,1)==0  % makes sense only for NON-DUPLICATES
                % Find the cell
                row=pc_poscell(ii,1);
                column=pc_poscell(ii,2);
                % The number of points in cell
                tc=mat_points_in_cell(row, column);
                % Generate vector of point IDs
                vect_aux=zeros(tc,1);
                vect_aux(1:tc)=grid_pts(row, column, 1:tc);
                % Take all the points corresponding to the cell
                for xx=1:tc  
                    jj=vect_aux(xx,1);  % Find the ID of the point
                    if jj>ii % makes sense to search duplicates only for GREATER IDs
                        if pc_in_xy (ii,1)==pc_in_xy (jj,1) && pc_in_xy (ii,2)==pc_in_xy(jj,2)
                            mat_control (jj,1)=1; % each duplicate gets an 1
                            index=index+1;
                        end  
                    end
                end
            end  
        end

        %% generate the output matrix (without duplicates)

        index_2=sum(mat_control); % count the number of duplicates
        nr_out=nr_vf-index_2;  % nr of single occurrences
        pc_out=zeros(nr_out,3);
        ind=0;

        for ii=1:nr_vf
            if mat_control (ii,1)==0  %% not a duplicate
               ind=ind+1; % a new 'true' record was found
               pc_out(ind,1:3)=pc_in(ii,1:3);  % copy the record on the appropriate row
            end
        end


        %disp(ind);
        end



