function [ Wgips , Cgips ] = Hexa_GlobalIPWeightsandCoords( RVEMODEL , Options )


%% PARSE GLOBAL DATA & FIRST COMPUTATIONS ======================================== %%

%-Get mesh data:
n_elem  = size(RVEMODEL.Conectivity,1);          %Number of Elements
npe     = size(RVEMODEL.Conectivity(:,3:end),2); %Nodes per element

%-Get quadrature data:
[ n_lips , Wlips , Clips ] = QuadratureData3D( Options.quadrature );

%-Compute local interpolating data:
switch npe
    case 8
        %elem_type = "3DHexa8N";
        InterpData = LagrangeInterp3D_C08N_Local(Clips);
    case 27
        %elem_type = "3DHexa27N";
        InterpData = LagrangeInterp3D_C127N_Local(Clips);
end

%-Global IP Conectivity : (WE ASSUME ALL ELEMENTS USE THE SAME QUADRATURE)
GIPsCon = [repelem((1:n_elem)',n_lips,1),repmat((1:n_lips)',n_elem,1)]; %[ Element N° , Local IP N° ]

%-Compute Interpolation at all IPs:
[ N_gips , detJ_gips , ~ , ~ , ~ ] = LagrangeInterp3D_Global_V( GIPsCon , InterpData , RVEMODEL );
if any(detJ_gips<=0.0)
    warning('Found non-positive Jacobian determinant');
    disp(detJ_gips(detJ_gips<=0.0));
end

%-Compute global weights for all ips:
Wgips = Wlips(GIPsCon(:,2))'.*detJ_gips;

%Compute cartesian coordinates of ips according to isoparametric interpolation: 
ncoords_x = RVEMODEL.Coordinates(:,1+1);
ncoords_y = RVEMODEL.Coordinates(:,1+2);
ncoords_z = RVEMODEL.Coordinates(:,1+3);
ncon_gips = permute( RVEMODEL.Conectivity(GIPsCon(:,1),3:end) , [1,3,2] ); 
            
Cgips = sum( N_gips.*[ncoords_x(ncon_gips) ,...
                      ncoords_y(ncon_gips) ,...
                      ncoords_z(ncon_gips) ] , 3 ); %[ X , Y , Z ]

end



