function [  ] = generateGraphicsCurvature (strct_elev, strct_gc, strct_mc, typ_curvature, met_curb, level_resol, bool_noise, level_noise, az, elev, shading_style, draw_colorbar, draw_figure)

%% Description
% Provides the graphical representation for a single curvature. The
% elevations are provided by the regular grid matrix, while the colors are
% provided by the associated curvatures.
% Input.
%   mat_elev: elevations of the regular grid
%   strct_gc, strct_mc: structures with values of curvatures for different levels of resolution
%   type_curvature: decide whether GC or MC
%   met_curb: method of curvature
%   az, elev: azimuth, elevation - parameters for 3D view
% Output. The grahical representation

hold off

%% Extract data from structures
% Curvature (GC or MC)
if typ_curvature==2
    strct=strct_mc;
else
    strct=strct_gc;
end
%% If no noise - the structure has the information regarding the levels of resolution
if bool_noise==0
    % Level of resolution as char
    level_char=['level' num2str(level_resol)];
    % Extract from the structure the grid / suitable column
    mat_elev=strct_elev.(matlab.lang.makeValidName(level_char));
    vect_col=strct.(matlab.lang.makeValidName(level_char))(:,met_curb);    
end

%% If noise - the structure has the information regarding the levels of noise
if bool_noise==1
    noise_char=['levelNoise' num2str(level_noise)];
    % Extract from the structure the grid / suitable column
    mat_elev=strct_elev.(matlab.lang.makeValidName(noise_char));
    vect_col=strct.(matlab.lang.makeValidName(noise_char))(:,met_curb);
end

%% Initialization
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
if draw_figure==1
    figure
end
surf (mat_elev,mat_col);
switch shading_style
    case 1 
        shading faceted
    case 2
        shading flat
    case 3
        shading interp
end
if draw_colorbar==0
   colorbar('off');
else
    colorbar;
end
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
        disp('Error!!!')
    end
    index=0;
    for ii=1:nr
        for jj=1:nc
         index=index+1;
         mat_out(ii,jj)=vect_in(index);
        end
    end

    end

