%% SCRIPT_3: 
%  Handle vertices and the computation of GC and MC (Dyn; Meyer) 



%% Initialize
% Angular defect
V_AD=zeros(nr_vf,2);
% GAUSSIAN CURVATURE
V_GC_1=zeros(nr_vf, 1); % GC 1: Dyn et al 2002
V_GC_2=zeros(nr_vf, 1); % GC 2: Meyer et al 2002
% MEAN CURVATURE
V_MC_1=zeros(nr_vf,1); % MC 1: Dyn et al 2002
V_MC_2=zeros(nr_vf,1); % MC 2: Meyer et al 2002

%% Main loop

for i=1:nr_vf
    %ANGLE DEFECT
    V_AD(i,1)=2*pi-V_AD_aux(i,1);
    V_AD(i,2)=(V_AD(i,1)*180)/pi;
    % GAUSSIAN CURVATURES 
    V_GC_1(i,1)=(3*V_AD(i,1))/V_area(i,1);
    V_GC_2(i,1)=V_AD(i,1)/V_mixedarea(i,1);
    % MEAN CURVATURES
    V_MC_1 (i,1)= (3*V_edge_contrib (i,1))/ V_area(i,1);
    V_MC_2 (i,1)= norm(V_mean_curvature_vector_aux(i,:))/(4*V_mixedarea(i,1));
end
 