function [ N , dNdX , detJ , J , invJ ] = LagrangeInterp3D_C127N( n_coords , nat_coords )
    
% ----------------------------------------------------------------------- %
% This function evaluates, for a given Gauss Point, the C1 sh.functions,  %
% its cartesian derivatives and the Jacobian Determinant for a 3D         %
% quadrature.                                                             %
%                                                                         %
% INPUTS:                                                                 %
%   [nat_coords] Local Coords of the point to evaluate the magnitudes     %
%   [n_coords]   Global Nodal Coordinates for the element                 %
%                                                                         %
% OUTPUTS:                                                                %
%   [N]     Vector with the Sh. functions evaluated at the given GP       %
%   [dNdX]  Vector w/the cartesian derivatives evaluated at the given GP  %
%   [detJ]  Jacobian Determinant evaluated at the given GP                %
%                                                                         %
% Additional Help - I/O Arrangement:                                      %
%   [ip_coords] = [ Eps_ip ; Eta_ip ; Dseta_ip ];                         %
%   [n_coords] = [ X1 Y1 Z1 ; .. .. .. ; X27 Y27 Z27 ];                   %
%   [ N ]  = [ N1 ; N2 ; ... ; N27 ];                                     %
%   [dNdX] = [ dN1/dX , dN1/dY , dN1/dZ;                                  %
%                    ...        ;                                         %
%              dN27/dX , dN27/dY , dN27/dZ ];                             %
%                                                                         %
% ------------------------------------------------- PWierna 04.V.2023 --- %


%% PARSE INPUT =================================================================== %%

%% Isoparametric nodal coordinates

pos_nodal = [-1 -1 -1; 1 -1 -1; 1  1 -1;-1  1 -1;...
             -1 -1  1; 1 -1  1; 1  1  1;-1  1  1;...
              0 -1 -1; 1  0 -1; 0  1 -1;-1  0 -1;...
             -1 -1  0; 1 -1  0; 1  1  0;-1  1  0;...
              0 -1  1; 1  0  1; 0  1  1;-1  0  1;...
              0 0 -1;...
              0 -1  0; 1  0  0; 0  1  0;-1  0  0;...
              0 0 1;...
              0 0 0];
t3 = nat_coords(:); %[ndim,1]
t4 = pos_nodal'; %[ndim,nodes]
t1 = 1 - t3.^2; %[ndim,1]
t2 = t4+t3; %[ndim,nodes]

v1 = 0.5.*t2.*t3.*(t4~=0) + t1.*(t4==0);
v2 = 0.5.*(t2+t3).*(t4~=0) - 2*t3.*(t4==0);

%% EVALUATE SHAPE FUNCTIONS ====================================================== %%
% t3 = nat_coords(:)';
% t4 = pos_nodal;
% t1 = 1 - t3.^2;
% t2 = t4+t3;
% N = prod(0.5.*t2.*t3.*(t4~=0) + t1.*(t4==0),2);

N = prod(v1,1)';

%% EVALUATE SHAPE FUNCTIONS LOCAL DERIVATIVES ==================================== %%

dN = v1([2,3,1],:).*v1([3,1,2],:).*v2;
dNdXi = dN(1,:)';
dNdEta = dN(2,:)';
dNddseta = dN(3,:)';

%% COMPUTE JACOBIAN, JAC DETERMINANT AND INVERSE ================================= %%
J = [ dNdXi , dNdEta , dNddseta ]' * n_coords;

detJ = J(1,1)*J(2,2)*J(3,3) - J(1,1)*J(2,3)*J(3,2) - J(1,2)*J(2,1)*J(3,3) +...
       J(1,2)*J(2,3)*J(3,1) + J(1,3)*J(2,1)*J(3,2) - J(1,3)*J(2,2)*J(3,1);      %Output

invJ = (1/detJ) * [ J(2,2)*J(3,3) - J(2,3)*J(3,2) , J(1,3)*J(3,2) - J(1,2)*J(3,3) , J(1,2)*J(2,3) - J(1,3)*J(2,2) ;
                    J(2,3)*J(3,1) - J(2,1)*J(3,3) , J(1,1)*J(3,3) - J(1,3)*J(3,1) , J(1,3)*J(2,1) - J(1,1)*J(2,3) ;
                    J(2,1)*J(3,2) - J(2,2)*J(3,1) , J(1,2)*J(3,1) - J(1,1)*J(3,2) , J(1,1)*J(2,2) - J(1,2)*J(2,1) ];


%% COMPUTE CARTESIAN DERIVATIVES OF THE SHAPE FUNCTIONS ========================== %%

% dNdX = [ dN1/dX , dN1/dY , dN1/dZ ; ... ; dN8/dX , dN8/dY , dN8/dZ ];
dNdX = ( invJ * [ dNdXi , dNdEta , dNddseta ]' )'; 

end