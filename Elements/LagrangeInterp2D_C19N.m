function [ N , dNdX , detJ , J , invJ ] = LagrangeInterp2D_C19N( n_coords , nat_coords )
    
% ----------------------------------------------------------------------- %
% This function evaluates, for a given Gauss Point, the C1 sh.functions,  %
% its cartesian derivatives and the Jacobian Determinant for a 2D         %
% quadrature.                                                             %
%                                                                         %
% INPUTS:                                                                 %
%   [nat_coords] Local Coords of the point to evaluate the magnitudes     %
%   [n_coords]   Global Nodal Coordinates for the element                 %
%                                                                         %
% OUTPUTS:                                                                %
%   [N]     Vector with the Sh. functions evaluated at the given GP       %
%   [dNdX]  Matrix w/the cartesian derivatives evaluated at the given GP  %
%   [detJ]  Jacobian Determinant evaluated at the given GP                %
%   [J]     Jacobian matrix evaluated at the given GP                     %
%                                                                         %
% Additional Help - I/O Arrangement:                                      %
%   [ip_coords] = [ Eps_ip ; Eta_ip ];                                    %
%   [n_coords] = [ X1 Y1 ; .. .. ; X9 Y9 ];                               %
%   [ N ]  = [ N1 ; N2 ; N3 ; ... ; N9 ];                                 %
%   [dNdX] = [ dN1/dX , dN1/dY ; ... ; dN9/dX , dN9/dY];                  %
%                                                                         %
% ----------------------------------------------- PWierna 29.VII.2022 --- %


%% PARSE INPUT ========================================================= %%

%Natural coordinates of the GP:
Xi = nat_coords(1);
Eta = nat_coords(2);


%% EVALUATE SHAPE FUNCTIONS ============================================ %%
N = zeros(9,1);

N(1) = (1/4) * (Xi^2 - Xi) * (Eta^2 - Eta);
N(3) = (1/4) * (Xi^2 + Xi) * (Eta^2 - Eta);
N(5) = (1/4) * (Xi^2 + Xi) * (Eta^2 + Eta);
N(7) = (1/4) * (Xi^2 - Xi) * (Eta^2 + Eta);

N(6) = (1/2) * (Eta^2 + Eta) * (1 - Xi^2); 
N(8) = (1/2) * (Xi^2 - Xi) * (1 - Eta^2); 
N(2) = (1/2) * (Eta^2 - Eta) * (1 - Xi^2); 
N(4) = (1/2) * (Xi^2 + Xi) * (1 - Eta^2); 

N(9) = (1 - Xi^2) * (1 - Eta^2); 

%% EVALUATE LOCAL DERIVATIVES ========================================== %%
DXi = zeros(9,1); %Vector with derivatives wrt Chi

DXi(1) = (1/4) * (Eta^2 - Eta)*(2*Xi - 1);
DXi(3) = (1/4) * (Eta^2 - Eta)*(2*Xi + 1);
DXi(5) = (1/4) * (Eta^2 + Eta)*(2*Xi + 1);
DXi(7) = (1/4) * (Eta^2 + Eta)*(2*Xi - 1);

DXi(6) = - (Eta^2 + Eta) * Xi;
DXi(8) = (1/2) * (1 - Eta^2)*(2*Xi - 1);
DXi(2) = - (Eta^2 - Eta) * Xi;
DXi(4) = (1/2) * (1 - Eta^2)*(2*Xi + 1);

DXi(9) = - 2*Xi * (1 - Eta^2);



DEta = zeros(9,1); %Vector with derivatives wrt Eta

DEta(1) = (1/4) * (Xi^2 - Xi)*(2*Eta - 1);
DEta(3) = (1/4) * (Xi^2 + Xi)*(2*Eta - 1);
DEta(5) = (1/4) * (Xi^2 + Xi)*(2*Eta + 1);
DEta(7) = (1/4) * (Xi^2 - Xi)*(2*Eta + 1);

DEta(6) = (1/2) * (1 - Xi^2)*(2*Eta + 1);
DEta(8) = - (Xi^2 - Xi) * Eta;
DEta(2) = (1/2) * (1 - Xi^2)*(2*Eta - 1);
DEta(4) = - (Xi^2 + Xi) * Eta;

DEta(9) = - 2*Eta * (1 - Xi^2);


%% EVALUATE ELEMENTS OF THE JACOBIAN =================================== %%
J = zeros(2);

J(1,1) = DXi' * n_coords(:,1);
J(1,2) = DXi' * n_coords(:,2);
J(2,1) = DEta' * n_coords(:,1);
J(2,2) = DEta' * n_coords(:,2);


%% COMPUTE JACOBIAN DETERMINANT AND INVERSE ============================ %%
detJ = J(1,1)*J(2,2) - J(2,1)*J(1,2); %Output
invJ = (1/detJ) * [J(2,2) -J(1,2); -J(2,1) J(1,1)];


%% COMPUTE CARTESIAN DERIVATIVES OF THE SHAPE FUNCTIONS ================ %%
%Assemble Transformation Matrix:
dNdX = zeros( 9 , 2 ); %=[ dN1/dX , dN1/dY ; ... ; dN9/dX , dN9/dY];

dNdX(1,:) = invJ * [ DXi(1) ; DEta(1) ];
dNdX(2,:) = invJ * [ DXi(2) ; DEta(2) ];
dNdX(3,:) = invJ * [ DXi(3) ; DEta(3) ];
dNdX(4,:) = invJ * [ DXi(4) ; DEta(4) ];
dNdX(5,:) = invJ * [ DXi(5) ; DEta(5) ];
dNdX(6,:) = invJ * [ DXi(6) ; DEta(6) ];
dNdX(7,:) = invJ * [ DXi(7) ; DEta(7) ];
dNdX(8,:) = invJ * [ DXi(8) ; DEta(8) ];
dNdX(9,:) = invJ * [ DXi(9) ; DEta(9) ];


return