function [ ] = generateGraphicsCorrelations(mat_data, vect_resol, mat_methods, vect_alphas)

%% Description
% Provides the graphical representation for the relationship between metrics,
% as quantified by the correlation coefficients / absolute errors.  
% Input.
%   mat_data: matrix to be represented (correlation coefficients / absolute
%   error)
%   vect_resol: vector of resolutions
%   mat_methods: information regarding the methods to be considered / the
%   reference method
%   vect_alphas: the vector of alpha values
% Output. The grahical representation
 
 

%% Initializations
[~,nr_alphas]=size(vect_alphas);
    aux=sum(vect_alphas);
    
%% Matrix of colors 
mat_colors=...
    [0.6 0.0 0.0;
    0.9 0.1 0.0; 
    1.0 0.6 0.;  
    0.8 0.6 0.0;
    0.0 0.8 0.0;
    0.1 0.9 1.0;
	0.1 0.5 0.4];
for ii=1:nr_alphas
    greyScale=vect_alphas(1,ii)/aux;
    mat_colors(7+ii,:)=[greyScale, greyScale,greyScale];
end
    mat_colors(8+nr_alphas,:)=[0 0 0];
 
%% Matrix of legends
mat_legend(1,:)='GB1     ';
mat_legend(2,:)='GB2     ';
mat_legend(3,:)='ET      ';
mat_legend(4,:)='TA      ';
mat_legend(5,:)='JF      ';
mat_legend(6,:)='NC(1R)  ';
mat_legend(7,:)='NC(2R)  '; 
for ii=1:nr_alphas
    aux=vect_alphas(1,ii);
    mat_legend(7+ii,:)=['IA(' sprintf('%.2f', aux) ')'];
end
mat_legend(8+nr_alphas,:)='RG/S    ';
   
 
%% Exclude the values corresponding to the reference method 
[~,nrMethods]=size(mat_methods);
for ii=1:nrMethods
    if mat_methods(2,ii)==1
        metRef=ii;
        referenceMethod=mat_legend(ii,:);
    end  
end

%% Exclude the colors and the legend items corresponding to methods not used
for ii=nrMethods:-1:1
    if mat_methods(1,ii)==0 || mat_methods(2,ii)==1
        mat_colors(ii,:)=[];
        mat_legend(ii,:)=[];
    end
end
sum_aux=sum(mat_methods(1,1:metRef-1));
columnRef=sum_aux+1;
mat_data(:,columnRef)=[];

%% Plot the data
hold off
[nr_levels_resolution, nr_represented_methods]=size(mat_data);
x=transpose(1:nr_levels_resolution);
y=mat_data;
xq=transpose(0.5:0.1:nr_levels_resolution+0.5);
yq=interp1(x,y,xq,'linear');
p1=plot(x,y,'o', 'LineWidth', 2);
for ii=1:nr_represented_methods
  p1(ii).Color=mat_colors(ii,:);  
end
hold on
p2=plot(xq,yq, 'LineWidth', 1.5);
for ii=1:nr_represented_methods
  p2(ii).Color=mat_colors(ii,:);  
end
warning off
set(gca,'XTick',1:nr_represented_methods,'XTickLabels', vect_resol);
leg=legend(mat_legend);
xlabel('Level of resolution');

end

