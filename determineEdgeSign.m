function [signEdge] = determineEdgeSign(P, Q, R1, R2)

%% Description
% Establishes the convexity of an edge for TERRAIN DATA
% Input. 
%   P, Q: determine the edge 
%   R1, R2: determine the adjacent faces
% Output.
%   signEdge: 1 if the edge is convex (0, 0, -1) inside or on the
%       polyhedral cone determined by the faces; -1 is the edge is concave; 0 if
%       the points are in the same plane
 

%% Initialization 
signEdge=1;

%% Coefficients of the plane equation
[~, ~, C1, ~]=determinePlaneEquation(P, Q, R1);
[~, ~, C2, ~]=determinePlaneEquation(P, Q, R2);

%% Signs of the complementary point / of the vector (0, 0, -1)
sign2=signHyperplane(P, Q, R1, R2); % position of R2 wrt (P Q R1)
sign2Vert=-C1; % position of (0, 0, -1) wrt the same plane (P Q R1)

sign1=signHyperplane(P, Q, R2, R1); % pos of R1 wrt (P Q R2)
sign1Vert=-C2; % position of (0, 0, -1) wrt the same plane (P Q R2)

%% Main test 
if sign(sign2)~=sign(sign2Vert) || sign(sign1)~=sign(sign1Vert)
    signEdge=-1;
end

%% Handle exceptions
if sign1==0 || sign2==0
    signEdge=0; % for coplanar points get 0   
end

if sign1Vert==0 || sign2Vert==0
    signEdge=0; % for vertical planes, by convention, the sign is 0
end

 

end

%% FUNCTIONS USED
    function [A, B, C, D] = determinePlaneEquation(P0, P1, P2)
    %% Determines the implicit equation of the plane passing through P0 P1 P2
    % Input: Points P0, P1, P2, represented as row vectors
    % Output: the coefficients (A B C D) of the plane equation Ax+By+Cz+D=0

    P0P1=P1-P0;
    P0P2=P2-P0;

    v = cross (P0P1, P0P2);
    A=v(1,1);
    B=v(1,2);
    C=v(1,3);
    D=-A*P0(1,1)-B*P0(1,2)-C*P0(1,3);
    end
    
    
    function [signHyp] = signHyperplane(P0, P1, P2, P3)

    %% Description
    % Determines the sign of a point P3 when replacing its coordinates in
    % the plane determined by the points P0, P1, P2

    [A, B, C, D]=determinePlaneEquation(P0, P1, P2);
    aux=A*P3(1,1)+B*P3(1,2)+C*P3(1,3)+D;
    signHyp=sign(aux);

    end

