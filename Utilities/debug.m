%-Get mesh data:
n_elem  = size(RVEMODEL.Conectivity,1);          %Number of Elements
npe     = size(RVEMODEL.Conectivity(:,3:end),2); %Nodes per element


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
n_gips  = n_lips * n_elem; %Total number of global ips 
GIPsCon = [repelem((1:n_elem)',n_lips,1),repmat((1:n_lips)',n_elem,1)]; %[ Element N° , Local IP N° ]

%-Compute Interpolation at all IPs:
[ ~ , detJ_gips , ~ , ~ , dNdX_gips ] = LagrangeInterp3D_Global_V( GIPsCon , InterpData , RVEMODEL );



all(all(abs(data_displ.strains - data_str.strains )<=1e-8))
all(all(abs(data_displ.stresses - data_str.stresses )<=1e-7))
all(all(abs(data_displ.fluctstrains - data_str.fluctstrains )<=1e-8))


fstrd=ComputeCurrentStrains_at_ips( data_displ.displacements , dNdX_gips , GIPsCon , RVEMODEL );
fstrd2=ComputeCurrentStrains_at_ips( d_d , dNdX_gips , GIPsCon , RVEMODEL );
fstrs=ComputeCurrentStrains_at_ips( data_str.displacements , dNdX_gips , GIPsCon , RVEMODEL );
fstrs2=ComputeCurrentStrains_at_ips( d_s , dNdX_gips , GIPsCon , RVEMODEL );

all(all(abs(fstrd - fstrd2)<=1e-8))
all(all(abs(fstrs - fstrs2)<=1e-8))
all(all(abs(fstrs - fstrd)<=1e-8))
all(all(abs(fstrs2 - fstrd2)<=1e-8))

mstrd = ComputeCurrentStrains_at_ips( data_displ.displmacro , dNdX_gips , GIPsCon , RVEMODEL );


function Strains_gips = ComputeCurrentStrains_at_ips( U , dNdX_gips , GIPsCon , MODEL )

%Hard setting:
dofpn = 3;

%Parse displacements:
Ux = U(1:dofpn:end);
Uy = U(2:dofpn:end);
Uz = U(3:dofpn:end);

%Parse shape functions derivatives:
dNdx_gips = permute( dNdX_gips(:,1,:) , [1,3,2] );
dNdy_gips = permute( dNdX_gips(:,2,:) , [1,3,2] );
dNdz_gips = permute( dNdX_gips(:,3,:) , [1,3,2] );

%Get global Node numbers for each ip:
ncon_gips = MODEL.Conectivity( GIPsCon(:,1) , 3:end); 

%Compute Voigt Strain vectors for all ips:
Strains_gips = zeros( size(GIPsCon,1) , 6 );
Strains_gips(:,1:3) = [ sum(dNdx_gips.*Ux(ncon_gips),2) , sum(dNdy_gips.*Uy(ncon_gips),2) , sum(dNdz_gips.*Uz(ncon_gips),2) ]; %[xx,yy,zz]
Strains_gips(:,4)   = sum(dNdy_gips.*Ux(ncon_gips),2) + sum(dNdx_gips.*Uy(ncon_gips),2); %[xy] 
Strains_gips(:,5)   = sum(dNdz_gips.*Ux(ncon_gips),2) + sum(dNdx_gips.*Uz(ncon_gips),2); %[xz] 
Strains_gips(:,6)   = sum(dNdz_gips.*Uy(ncon_gips),2) + sum(dNdy_gips.*Uz(ncon_gips),2); %[yz] 

end