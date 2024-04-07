function [ InterpData , FacesData ] = LagrangeInterp3D_C127N_Local( Natural_Coords )

%Get number of integration points:
n_ips = size(Natural_Coords,1);

%Initialize array:
InterpData = zeros(27,4,n_ips);

   %Note: Output Data Arrangement
   % Rows = Shape F. Number (Node)
   % Cols = (:,1,:): Sh.Functions [N]
   %        (:,2,:): Derivatives wrt Chi [dNdXi]
   %        (:,3,:): Derivatives wrt Eta [dNdEta]
   %        (:,4,:): Derivatives wrt Dseta [dNddseta]
   % 3rd Dim = Local IP Number (ip)


pos_nodal = [-1 -1 -1; 1 -1 -1; 1  1 -1;-1  1 -1;...
             -1 -1  1; 1 -1  1; 1  1  1;-1  1  1;...
              0 -1 -1; 1  0 -1; 0  1 -1;-1  0 -1;...
             -1 -1  0; 1 -1  0; 1  1  0;-1  1  0;...
              0 -1  1; 1  0  1; 0  1  1;-1  0  1;...
              0 0 -1;...
              0 -1  0; 1  0  0; 0  1  0;-1  0  0;...
              0 0 1;...
              0 0 0];
          
t3 = permute(Natural_Coords,[3,2,1]); %[1,dim,gp]
t4 = pos_nodal; %[nodes,dim,1]
t1 = 1 - t3.^2; %[1,dim,gp]
t2 = t4+t3; %[nodes,dim,gp]

v1 = 0.5.*t2.*t3.*(t4~=0) + t1.*(t4==0); %[nodes,dim,gp]
v2 = 0.5.*(t2+t3).*(t4~=0) - 2*t3.*(t4==0); %[nodes,dim,gp]
   
InterpData(:,1,:) = prod(v1,2);
InterpData(:,2:4,:) = v1(:,[2,3,1],:).*v1(:,[3,1,2],:).*v2;

%Faces Data: Local Faces Normals & Connectivities (Row = Local Face number)
FacesData.Normals = [ -1,0,0 ; 1,0,0 ; 0,-1,0 ; 0,1,0 ; 0,0,-1 ; 0,0,1 ];
FacesData.Connectivities = [1 4 5 8 12 13 16 20 25;
                            2 3 6 7 10 14 15 18 23;
                            1 2 5 6 9 13 14 17 22;
                            3 4 7 8 11 15 16 19 24;
                            1 2 3 4 9 10 11 12 21;
                            5 6 7 8 17 18 19 20 26 ];
                                         
end